importall LegHemodynamics

function main()

rstflag = "yes" # restarting from scratch or previous simulation
hemoflag = "no" # 10% hemorrhage from left femoral artery
saveflag = "yes" # save solution to .mat file
coupleflag = "no" # coupling to 3D liver tissue model
timeflag = "yes" # solver timing
assimflag = "yes" # patient data assimilation via EnKF

# build solution struct or generate ensemble
if assimflag == "no"
    if rstflag == "no"
        loadfile = "arterylist_2.txt"; # default artery data file for new sim
    elseif rstflag == "yes"
        loadfile = "test.mat"; # restart file
    end
    system = LegHemodynamics.buildall(loadfile;numbeatstotal=1,restart=rstflag);
    savefile = "test.mat" # filename for saving (only used if saveflag == "yes")
elseif assimflag == "yes"
    ensemblesize = 10;
    if rstflag == "no"
        loadfiles = ["arterylist_2_thick_branches.txt" for i=1:ensemblesize];
    elseif rstflag == "yes"
        loadfiles = ["small_dev_7_$i.mat" for i=1:ensemblesize];
    end
    systems = pmap((a1)->LegHemodynamics.buildall(a1;numbeatstotal=1,restart=rstflag),loadfiles);
    savefiles = ["converged_small_dev_1_$i.mat" for i=1:ensemblesize];
    if rstflag == "no" # randomize ICs on new ensemble
        soln = pmap((a1)->LegHemodynamics.applycustomics!(a1),systems);
        systems = [soln[i] for i=1:ensemblesize];
    elseif rstflag == "yes" # ensure uniform time step
        h = [];
        for i = 1:length(systems)
            push!(h,systems[i].solverparams.h)
        end
        for i = 1:length(systems)
            systems[i].solverparams.h = minimum(h);
        end
        systems = pmap((a1,a2)->LegHemodynamics.rediscretizet!(a1,a2),systems,zeros(Int64,length(systems)));
    end
end

# create interpolation object for upstream flow velocity data (needed for proximal BC)
patientdata = MAT.matread("flowrates_mL.mat");
tdata = reshape(patientdata["t"],101);
udata = reshape(patientdata["u_pop_cms"],101)*LegHemodynamics.cmTom;
if assimflag == "no"
    spl = Dierckx.Spline1D(tdata,udata;k=1);
elseif assimflag == "yes"
    spl = [Dierckx.Spline1D(tdata,udata;k=1) for i=1:ensemblesize];
    udata = reshape(patientdata["u_ant_cms"],101)*LegHemodynamics.cmTom;
    splant = Dierckx.Spline1D(tdata,udata;k=1);
    udata = reshape(patientdata["u_post_cms"],101)*LegHemodynamics.cmTom;
    splpost = Dierckx.Spline1D(tdata,udata;k=1);
    udata = reshape(patientdata["u_per_cms"],101)*LegHemodynamics.cmTom;
    splper = Dierckx.Spline1D(tdata,udata;k=1);
end

# timers
if assimflag == "no"
    times = LegHemodynamics.CVTimer();
elseif assimflag == "yes"
    times = [LegHemodynamics.CVTimer() for i=1:ensemblesize];
end

# collect all IDs of terminal/non-terminal branches
terms = Int64[];
splits = Int64[];
if assimflag == "no"
    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            push!(terms,i);
        else
            push!(splits,i);
        end
    end
elseif assimflag == "yes"
    for i = 1:length(systems[1].branches.ID)
        if isempty(systems[1].branches.children[i])
            push!(terms,i);
        else
            push!(splits,i);
        end
    end
    term_itr = [terms for i=1:ensemblesize];
    split_itr = [splits for i=1:ensemblesize];
end

# ensemble distributions, allocators, and scalings
if assimflag == "yes"
    nparams = 15;
    nstates = 3;
    nmeas = 3;
    errors = LegHemodynamics.Errors(nparams);

    X = [zeros(nstates) for i=1:ensemblesize];
    θ = [zeros(nparams) for i=1:ensemblesize];
    Y = [zeros(nmeas) for i=1:ensemblesize];

    xhat = zeros(nstates);
    θhat = zeros(nparams);
    yhat = zeros(nmeas);

    yi = [zeros(nmeas) for i=1:ensemblesize];

    # parameter scalings
    θs = zeros(nparams);
    for i=1:ensemblesize
        θs[1] += systems[i].branches.term[14].R[2]; # anterior
        θs[2] += systems[i].branches.term[24].R[2]; # peroneal
        θs[3] += systems[i].branches.term[32].R[2]; # posterior
        θs[4] += systems[i].branches.term[14].C[1];
        θs[5] += systems[i].branches.term[24].C[1];
        θs[6] += systems[i].branches.term[32].C[1];
        θs[7] += systems[i].branches.termscalings[1];
        θs[8] += systems[i].branches.termscalings[2];
        θs[9] += systems[i].branches.termscalings[3];
        θs[10] += systems[i].branches.termscalings[4];
        θs[11] += systems[i].branches.termscalings[5];
        θs[12] += systems[i].branches.termscalings[6];
        # θs[9] += systems[i].branches.A0[14];
        # θs[10] += systems[i].branches.A0[24];
        # θs[11] += systems[i].branches.A0[31];
        θs[13] += systems[i].branches.beta[14];
        θs[14] += systems[i].branches.beta[24];
        θs[15] += systems[i].branches.beta[32];
        # θs[9] += systems[i].branches.YModinPa[14];
        # θs[10] += systems[i].branches.YModinPa[24];
        # θs[11] += systems[i].branches.YModinPa[31];
    end
    θs /= ensemblesize;

    # RTPS allocators
    p = 0.5; # RTPS relaxation amount, ∃ [0.5,0.95]
    c = zeros(ensemblesize);
    σxb = zeros(nstates);
    σxa = zeros(nstates);
    σtb = zeros(nparams);
    σta = ones(nparams);

    tout = Float64[];
    xout = Vector{Float64}[];
    xoutv = Vector{Float64}[];
    yout = Vector{Float64}[];
    youtv = Vector{Float64}[];
    innov = Vector{Float64}[];
    ioutv = Vector{Float64}[];
    θout = Vector{Float64}[];
    θoutv = Vector{Float64}[];
    Pθout = Vector{Float64}[];
    Pθoutv = Vector{Float64}[];
    Pxout = Vector{Float64}[];
    Pxoutv = Vector{Float64}[];
    lbx = Vector{Float64}[];
    ubx = Vector{Float64}[];
    lbxv = Vector{Float64}[];
    ubxv = Vector{Float64}[];

    # default parameter values for non-measured terminal constraints
    Rdefault = 0.75*[77.6909,141.4980,114.7901,102.6618,176.7149,257.4914,18.6383,
        107.4268,291.8191,191.9240,149.3493,35.9008,53.0083,150.1083,259.9696,
        142.5643,95.1046]*LegHemodynamics.gTokg/(LegHemodynamics.mmTom^4);
    Cdefault = [0.08972235,0.04926296,0.06072485,0.06789874,0.03944552,0.02707124,
        0.39707822,0.06488707,0.02388675,0.03631963,0.04667319,0.19569034,0.13150028,
        0.04643721,0.02681318,0.04830696,0.00723676]*(LegHemodynamics.mmTom^4)/LegHemodynamics.gTokg;
end

# time step counter
if assimflag == "no"
    n = system.solverparams.nstart;
elseif assimflag == "yes"
    n = [systems[1].solverparams.nstart for i=1:ensemblesize];
end

# solver loop
tic();
if assimflag == "no"
    while system.solverparams.numbeats < system.solverparams.numbeatstotal
        # time within heart cycle
        tp = system.t[n+1] - sum(system.solverparams.th*system.solverparams.numbeats);
        if mod(n,1000) == 0
            println("Reached time step $n. t = $tp s.")
        end
        # interpolate proximal flow velocity from patient data
        uprox = spl(tp);
        # TVD RK3 time integration
        LegHemodynamics.tvdrk3!(system,times,n,splits,terms,-uprox);
        tic();
        # arterial pressure update
        LegHemodynamics.arterialpressure!(system,n);
        # check for start of new cardiac cycle
        oldnumbeats = system.solverparams.numbeats;
        LegHemodynamics.setnumbeats!(system,n);
        if oldnumbeats < system.solverparams.numbeats
            # update time discretization
            if system.solverparams.numbeats < system.solverparams.numbeatstotal
                LegHemodynamics.updatediscretization!(system);
            end
        end
        n+=1;
        times.tr += toq();
    end
    times.tt += toc();
elseif assimflag == "yes"
    nsamp = [300 for i=1:ensemblesize];
    while systems[1].solverparams.numbeats < systems[1].solverparams.numbeatstotal
        # normal time integration between measurements
        rflag = "predict";
        soln = pmap((a1,a2,a3,a4,a5,a6,a7)->LegHemodynamics.advancetime!(a1,a2,a3,a4,a5,a6,a7;
            runtype=rflag),systems,spl,times,n,nsamp,term_itr,split_itr);
        systems = [soln[i][1] for i=1:ensemblesize];
        n = [soln[i][2] for i=1:ensemblesize];
        times = [soln[i][3] for i=1:ensemblesize];
        # data assimilation
        if assimflag == "yes" && mod((n[1]-1)/nsamp[1],1) == 0
            println("Current time step: $(n[1])")
            println("Current time: $(systems[1].t[n[1]])")

            # time within heart cycle
            tp = systems[1].t[n[1]] - sum(systems[1].solverparams.th*systems[1].solverparams.numbeats);
            # non-dimensional measurement (velocity)
            y = -[splant(tp),splpost(tp),splper(tp)]/LegHemodynamics.us;
            # y = -[splant(tp)]/LegHemodynamics.us;
            println("Scaled velocity measurements at t = $tp s: $y")

            # # generate measurement replicates
            # for i = 1:length(systems)
            #     yi[i][1] = y[1] + rand(Distributions.Normal(0,errors.odev[1]/LegHemodynamics.us));
            #     yi[i][2] = y[2] + rand(Distributions.Normal(0,errors.odev[2]/LegHemodynamics.us));
            #     yi[i][3] = y[3] + rand(Distributions.Normal(0,errors.odev[3]/LegHemodynamics.us));
            # end
            #
            # # forecast parameters and means
            # for i = 1:length(systems)
            #     θ[i][1] = systems[i].branches.term[14].R[2];
            #     θ[i][2] = systems[i].branches.term[24].R[2];
            #     θ[i][3] = systems[i].branches.term[32].R[2];
            #     θ[i][4] = systems[i].branches.term[14].C[1];
            #     θ[i][5] = systems[i].branches.term[24].C[1];
            #     θ[i][6] = systems[i].branches.term[32].C[1];
            #     θ[i][7] = systems[i].branches.termscalings[1];
            #     θ[i][8] = systems[i].branches.termscalings[2];
            #     θ[i][9] = systems[i].branches.termscalings[3];
            #     θ[i][10] = systems[i].branches.termscalings[4];
            #     θ[i][11] = systems[i].branches.termscalings[5];
            #     θ[i][12] = systems[i].branches.termscalings[6];
            #     # θ[i][9] = systems[i].branches.A0[14];
            #     # θ[i][10] = systems[i].branches.A0[24];
            #     # θ[i][11] = systems[i].branches.A0[31];
            #     θ[i][13] = systems[i].branches.beta[14];
            #     θ[i][14] = systems[i].branches.beta[24];
            #     θ[i][15] = systems[i].branches.beta[32];
            #     # θ[i][9] = systems[i].branches.YModinPa[14];
            #     # θ[i][10] = systems[i].branches.YModinPa[24];
            #     # θ[i][11] = systems[i].branches.YModinPa[31];
            # end
            #
            # # forecast parameters (dimensional)
            # θ = [LegHemodynamics.paramwalk!(errors,θs,θ[i]) for i in 1:length(systems)];
            #
            # # non-dimensionalize parameters
            # for i = 1:length(systems)
            #     θ[i] = θ[i]./θs;
            # end
            #
            # # RTPS prior standard deviation
            # for i in 1:length(θ[1])
            #     for j in 1:length(systems)
            #         c[j] = θ[j][i];
            #     end
            #     σtb[i] = std(c;corrected=true);
            # end
            #
            # # parameters back into model
            # for i = 1:length(soln)
            #     systems[i].branches.term[14].R[2] = θ[i][1]*θs[1];
            #     systems[i].branches.term[24].R[2] = θ[i][2]*θs[2];
            #     systems[i].branches.term[32].R[2] = θ[i][3]*θs[3];
            #     systems[i].branches.term[14].C[1] = θ[i][4]*θs[4];
            #     systems[i].branches.term[24].C[1] = θ[i][5]*θs[5];
            #     # systems[i].branches.term[24].C[1] = 2.5e-11;
            #     systems[i].branches.term[32].C[1] = θ[i][6]*θs[6];
            #     systems[i].branches.termscalings[1] = θ[i][7]*θs[7];
            #     systems[i].branches.termscalings[2] = θ[i][8]*θs[8];
            #     systems[i].branches.termscalings[3] = θ[i][9]*θs[9];
            #     systems[i].branches.termscalings[4] = θ[i][10]*θs[10];
            #     systems[i].branches.termscalings[5] = θ[i][11]*θs[11];
            #     systems[i].branches.termscalings[6] = θ[i][12]*θs[12];
            #     # k = 1;
            #     for j in [2,4;] # popliteal
            #         systems[i].branches.term[j].C[1] *= 1/3*(systems[i].branches.termscalings[1]+
            #             systems[i].branches.termscalings[3]+systems[i].branches.termscalings[5]);
            #         systems[i].branches.term[j].R[2] *= 1/3*(systems[i].branches.termscalings[2]+
            #             systems[i].branches.termscalings[4]+systems[i].branches.termscalings[6]);
            #         # systems[i].branches.term[j].C[1] = 1/3*(systems[i].branches.termscalings[1]+
            #             # systems[i].branches.termscalings[3]+systems[i].branches.termscalings[5])*Cdefault[k];
            #         # systems[i].branches.term[j].R[2] = 1/3*(systems[i].branches.termscalings[2]+
            #             # systems[i].branches.termscalings[4]+systems[i].branches.termscalings[6])*Rdefault[k];
            #         # k+=1;
            #     end
            #     for j in [7,9,11,13;] # anterior
            #         systems[i].branches.term[j].C[1] *= systems[i].branches.termscalings[1];
            #         systems[i].branches.term[j].R[2] *= systems[i].branches.termscalings[2];
            #         # systems[i].branches.term[j].C[1] = systems[i].branches.termscalings[1]*Cdefault[k];
            #         # systems[i].branches.term[j].R[2] = systems[i].branches.termscalings[2]*Rdefault[k];
            #         # k+=1;
            #     end
            #     for j in [17,19,21,23;] # posterior
            #         systems[i].branches.term[j].C[1] *= systems[i].branches.termscalings[3];
            #         systems[i].branches.term[j].R[2] *= systems[i].branches.termscalings[4];
            #         # systems[i].branches.term[j].C[1] = systems[i].branches.termscalings[3]*Cdefault[k];
            #         # systems[i].branches.term[j].R[2] = systems[i].branches.termscalings[4]*Rdefault[k];
            #         # k+=1;
            #     end
            #     for j in [25,27,29,31;] # peroneal
            #         systems[i].branches.term[j].C[1] *= systems[i].branches.termscalings[5];
            #         systems[i].branches.term[j].R[2] *= systems[i].branches.termscalings[6];
            #         # systems[i].branches.term[j].C[1] = systems[i].branches.termscalings[5]*Cdefault[k];
            #         # systems[i].branches.term[j].R[2] = systems[i].branches.termscalings[6]*Rdefault[k];
            #         # k+=1;
            #     end
            #     for j in [1:4;] # popliteal
            #         systems[i].branches.beta[j] *= 1/3*(θ[i][13]*θs[13]/systems[i].branches.beta[14]+
            #             θ[i][14]*θs[14]/systems[i].branches.beta[24]+θ[i][15]*θs[15]/systems[i].branches.beta[32]);
            #     end
            #     for j in [5,7:13;] # anteriortibial
            #         systems[i].branches.beta[j] *= θ[i][13]*θs[13]/systems[i].branches.beta[14];
            #     end
            #     for j in [6,15,17:23;] # posteriortibial
            #         systems[i].branches.beta[j] *= θ[i][14]*θs[14]/systems[i].branches.beta[24];
            #     end
            #     for j in [16,25:31;] # peroneal
            #         systems[i].branches.beta[j] *= θ[i][15]*θs[15]/systems[i].branches.beta[32];
            #     end
            #     systems[i].branches.beta[14] = θ[i][13]*θs[13];
            #     systems[i].branches.beta[24] = θ[i][14]*θs[14];
            #     systems[i].branches.beta[32] = θ[i][15]*θs[15];
            #     # systems[i].branches.A0[14] = θ[i][9]*θs[9];
            #     # systems[i].branches.A0[24] = θ[i][10]*θs[10];
            #     # systems[i].branches.A0[31] = θ[i][11]*θs[11];
            #     # for j in [14,24,32]
            #     #     systems[i].branches.term[j].R[1] = systems[i].solverparams.rho*
            #     #         systems[i].branches.c0[j][end]/systems[i].branches.A0[j];
            #     # end
            #     for j = 1:length(systems[i].branches.ID)
            #         systems[i].branches.c0[j][end] = sqrt(0.5*systems[i].branches.beta[j]/
            #             systems[i].solverparams.rho)*systems[i].branches.A0[j]^0.25
            #         if isempty(systems[i].branches.children[j])
            #             systems[i].branches.term[j].R[1] = systems[i].solverparams.rho*
            #                 systems[i].branches.c0[j][end]/systems[i].branches.A0[j];
            #         end
            #         # systems[i].branches.beta[j] = sqrt(pi)*systems[i].branches.thicknessinmm[j]*LegHemodynamics.mmTom*
            #             # systems[i].branches.YModinPa[j]/((1-0.5^2)*systems[i].branches.A0[j]);
            #     end
            # end
            #
            # # find new time step to satisfy CFL condition for all members
            # h = [];
            # for i = 1:length(soln)
            #     for j in 1:length(systems[i].branches.ID)
            #         push!(h,systems[i].solverparams.CFL*systems[i].branches.k[j]/
            #             systems[i].branches.c0[j][end]);
            #     end
            # end
            # if minimum(h) != systems[1].solverparams.h
            #     # println("Old time step size: $(systems[1].solverparams.h) s")
            #     for i = 1:length(soln)
            #         systems[i].solverparams.h = minimum(h);
            #     end
            #     # println("New time step size: $(systems[1].solverparams.h) s")
            #     # re-discretize
            #     if systems[1].solverparams.numbeats < systems[1].solverparams.numbeatstotal
            #         systems = pmap((a1,a2)->LegHemodynamics.rediscretizet!(a1,a2),systems,n);
            #     end
            # end
            #
            # # forecast step
            # rflag = "assim";
            # soln = pmap((a1,a2,a3,a4,a5,a6,a7)->LegHemodynamics.advancetime!(a1,a2,a3,a4,a5,a6,a7;
            #     runtype=rflag),systems,spl,times,n,nsamp,term_itr,split_itr);
            # systems = [soln[i][1] for i=1:ensemblesize];
            # n = [soln[i][2] for i=1:ensemblesize];
            # times = [soln[i][3] for i=1:ensemblesize];
            #
            # # forecast state, measurements
            # X = [zeros(nstates) for i=1:ensemblesize];
            # Y = [zeros(nmeas) for i=1:ensemblesize];
            #
            # for i = 1:ensemblesize
            #     X[i][1] = systems[i].branches.Q[14][1,n[1]]/
            #         (systems[i].branches.A[14][1,n[1]]*LegHemodynamics.us);
            #     X[i][2] = systems[i].branches.Q[24][1,n[1]]/
            #         (systems[i].branches.A[24][1,n[1]]*LegHemodynamics.us);
            #     X[i][3] = systems[i].branches.Q[32][1,n[1]]/
            #         (systems[i].branches.A[32][1,n[1]]*LegHemodynamics.us);
            #     Y[i][1] = systems[i].branches.Q[14][1,n[1]]/
            #         (systems[i].branches.A[14][1,n[1]]*LegHemodynamics.us);
            #     Y[i][2] = systems[i].branches.Q[24][1,n[1]]/
            #         (systems[i].branches.A[24][1,n[1]]*LegHemodynamics.us);
            #     Y[i][3] = systems[i].branches.Q[32][1,n[1]]/
            #         (systems[i].branches.A[32][1,n[1]]*LegHemodynamics.us);
            # end
            #
            # # forecast mean state, parameters, measurement
            # xhat = mean(X);
            # yhat = mean(Y);
            # θhat = mean(θ);
            # # println("Normalized ̂x, first forecast: $xhat")
            # println("Normalized ̂y, first forecast: $yhat")
            # println("Normalized ̂θ, first forecast: $θhat")
            #
            # # forecast params./meas. cross covariance, meas. covariance
            # Pty = zeros(length(θhat),length(yhat))
            # Pyy = zeros(length(yhat),length(yhat))
            # for i = 1:length(soln)
            #     Pty += *((θ[i] .- θhat),(Y[i] .- yhat)');
            #     Pyy += *((Y[i] .- yhat),(Y[i] .- yhat)');
            # end
            # Pty ./= (ensemblesize);
            # Pyy ./= (ensemblesize);
            #
            # # println("Parameter-measurement cross-covariance:")
            # # display(Pty)
            #
            # # add meas. noise to meas. covariance (allows invertibility)
            # Pyy += diagm([(errors.odev[1]/LegHemodynamics.us)^2,
            #     (errors.odev[2]/LegHemodynamics.us)^2,
            #     (errors.odev[3]/LegHemodynamics.us)^2],0);
            # # Pyy += diagm([(errors.odev[1]/LegHemodynamics.us)^2],0);
            #
            # # println("Measurement covariance:")
            # # display(Pyy)
            #
            # # parameter Kalman gain
            # K = Pty*inv(Pyy);
            # # println("Parameter Kalman gain:")
            # # display(K)
            #
            # # parameter analysis step
            # for i = 1:length(soln)
            #     # println("ith measurement replicate: $(yi[i])")
            #     # println("ith forecast measurement: $(Y[i])")
            #     θ[i][:] += K*(yi[i] .- Y[i]);
            # end
            #
            # # RTPS parameter covariance inflation
            # for i in 1:length(θ[1])
            #     for j in 1:length(systems)
            #         c[j] = θ[j][i];
            #     end
            #     σta[i] = std(c;corrected=true);
            # end
            # θhat = mean(θ);
            # for i = 1:length(soln)
            #     θ[i] .= θ[i] .+ p.*((σtb.-σta)./σta).*(θ[i].-θhat);
            # end
            #
            # # analysis parameters back into ensemble members
            # for i = 1:length(soln)
            #     for j = 1:nparams
            #         if θ[i][j]*θs[j] <= errors.lb[j]
            #             # println("Warning: analysis θ_$j below lower bound for member $i.
            #                 # Setting to lower bound of $(errors.lb[j]).")
            #             θ[i][j] = errors.lb[j]/θs[j];
            #         end
            #         if θ[i][j]*θs[j] >= errors.ub[j]
            #             # println("Warning: analysis θ_$j above upper bound for member $i.
            #                 # Setting to upper bound of $(errors.ub[j]).")
            #             θ[i][j] = errors.ub[j]/θs[j];
            #         end
            #     end
            #     systems[i].branches.term[14].R[2] = θ[i][1]*θs[1];
            #     systems[i].branches.term[24].R[2] = θ[i][2]*θs[2];
            #     systems[i].branches.term[32].R[2] = θ[i][3]*θs[3];
            #     systems[i].branches.term[14].C[1] = θ[i][4]*θs[4];
            #     systems[i].branches.term[24].C[1] = θ[i][5]*θs[5];
            #     # systems[i].branches.term[24].C[1] = 2.5e-11;
            #     systems[i].branches.term[32].C[1] = θ[i][6]*θs[6];
            #     systems[i].branches.termscalings[1] = θ[i][7]*θs[7];
            #     systems[i].branches.termscalings[2] = θ[i][8]*θs[8];
            #     systems[i].branches.termscalings[3] = θ[i][9]*θs[9];
            #     systems[i].branches.termscalings[4] = θ[i][10]*θs[10];
            #     systems[i].branches.termscalings[5] = θ[i][11]*θs[11];
            #     systems[i].branches.termscalings[6] = θ[i][12]*θs[12];
            #     # k = 1;
            #     for j in [2,4;] # popliteal
            #         systems[i].branches.term[j].C[1] *= 1/3*(systems[i].branches.termscalings[1]+
            #             systems[i].branches.termscalings[3]+systems[i].branches.termscalings[5]);
            #         systems[i].branches.term[j].R[2] *= 1/3*(systems[i].branches.termscalings[2]+
            #             systems[i].branches.termscalings[4]+systems[i].branches.termscalings[6]);
            #         # systems[i].branches.term[j].C[1] = 1/3*(systems[i].branches.termscalings[1]+
            #             # systems[i].branches.termscalings[3]+systems[i].branches.termscalings[5])*Cdefault[k];
            #         # systems[i].branches.term[j].R[2] = 1/3*(systems[i].branches.termscalings[2]+
            #             # systems[i].branches.termscalings[4]+systems[i].branches.termscalings[6])*Rdefault[k];
            #         # k+=1;
            #     end
            #     for j in [7,9,11,13;] # anterior
            #         systems[i].branches.term[j].C[1] *= systems[i].branches.termscalings[1];
            #         systems[i].branches.term[j].R[2] *= systems[i].branches.termscalings[2];
            #         # systems[i].branches.term[j].C[1] = systems[i].branches.termscalings[1]*Cdefault[k];
            #         # systems[i].branches.term[j].R[2] = systems[i].branches.termscalings[2]*Rdefault[k];
            #         # k+=1;
            #     end
            #     for j in [17,19,21,23;] # posterior
            #         systems[i].branches.term[j].C[1] *= systems[i].branches.termscalings[3];
            #         systems[i].branches.term[j].R[2] *= systems[i].branches.termscalings[4];
            #         # systems[i].branches.term[j].C[1] = systems[i].branches.termscalings[3]*Cdefault[k];
            #         # systems[i].branches.term[j].R[2] = systems[i].branches.termscalings[4]*Rdefault[k];
            #         # k+=1;
            #     end
            #     for j in [25,27,29,31;] # peroneal
            #         systems[i].branches.term[j].C[1] *= systems[i].branches.termscalings[5];
            #         systems[i].branches.term[j].R[2] *= systems[i].branches.termscalings[6];
            #         # systems[i].branches.term[j].C[1] = systems[i].branches.termscalings[5]*Cdefault[k];
            #         # systems[i].branches.term[j].R[2] = systems[i].branches.termscalings[6]*Rdefault[k];
            #         # k+=1;
            #     end
            #     for j = 1:length(term_itr[1])
            #         if term_itr[i][j] != 14 && term_itr[i][j] != 24 && term_itr[i][j] != 32
            #             # systems[i].branches.term[term_itr[i][j]].C[1] *= systems[i].branches.termscalings[1];
            #             # if systems[i].branches.term[term_itr[i][j]].C[1] < 1e-11
            #             #     systems[i].branches.term[term_itr[i][j]].C[1] = 1e-11;
            #             # end
            #             if systems[i].branches.term[term_itr[i][j]].C[1] < 0.5*Cdefault[j]
            #                 # println("Warning: Analysis terminal C < lower bound for member $i,
            #                     # artery $j. Setting to 1e-11.")
            #                 systems[i].branches.term[term_itr[i][j]].C[1] = 0.5*Cdefault[j];
            #             elseif systems[i].branches.term[term_itr[i][j]].C[1] < 1e-12
            #                 systems[i].branches.term[term_itr[i][j]].C[1] = 1e-12;
            #             elseif systems[i].branches.term[term_itr[i][j]].C[1] > 2*Cdefault[j];
            #                 systems[i].branches.term[term_itr[i][j]].C[1] = 2*Cdefault[j];
            #             end
            #             # systems[i].branches.term[term_itr[i][j]].R[2] *= systems[i].branches.termscalings[2];
            #             if systems[i].branches.term[term_itr[i][j]].R[2] < 0.5*Rdefault[j]
            #                 systems[i].branches.term[term_itr[i][j]].R[2] = 0.5*Rdefault[j];
            #             elseif systems[i].branches.term[term_itr[i][j]].R[2] > 2*Rdefault[j]
            #                 systems[i].branches.term[term_itr[i][j]].R[2] = 2*Rdefault[j];
            #             end
            #         end
            #     end
            #     for j in [1:4;] # popliteal
            #         systems[i].branches.beta[j] *= 1/3*(θ[i][13]*θs[13]/systems[i].branches.beta[14]+
            #             θ[i][14]*θs[14]/systems[i].branches.beta[24]+θ[i][15]*θs[15]/systems[i].branches.beta[32]);
            #     end
            #     for j in [5,7:13;] # anteriortibial
            #         systems[i].branches.beta[j] *= θ[i][13]*θs[13]/systems[i].branches.beta[14];
            #     end
            #     for j in [6,15,17:23;] # posteriortibial
            #         systems[i].branches.beta[j] *= θ[i][14]*θs[14]/systems[i].branches.beta[24];
            #     end
            #     for j in [16,25:31;] # peroneal
            #         systems[i].branches.beta[j] *= θ[i][15]*θs[15]/systems[i].branches.beta[32];
            #     end
            #     systems[i].branches.beta[14] = θ[i][13]*θs[13];
            #     systems[i].branches.beta[24] = θ[i][14]*θs[14];
            #     systems[i].branches.beta[32] = θ[i][15]*θs[15];
            #     # systems[i].branches.A0[14] = θ[i][9]*θs[9];
            #     # systems[i].branches.A0[24] = θ[i][10]*θs[10];
            #     # systems[i].branches.A0[31] = θ[i][11]*θs[11];
            #     for j = 1:length(systems[i].branches.ID)
            #         # systems[i].branches.beta[j] = sqrt(pi)*systems[i].branches.thicknessinmm[j]*LegHemodynamics.mmTom*
            #             # systems[i].branches.YModinPa[j]/((1-0.5^2)*systems[i].branches.A0[j]);
            #         systems[i].branches.c0[j][end] = sqrt(0.5*systems[i].branches.beta[j]/
            #             systems[i].solverparams.rho)*systems[i].branches.A0[j]^0.25
            #         if isempty(systems[i].branches.children[j])
            #             systems[i].branches.term[j].R[1] = systems[i].solverparams.rho*
            #                     systems[i].branches.c0[j][end]/systems[i].branches.A0[j];
            #         end
            #     end
            #     # for j in [14,24,32]
            #     #     systems[i].branches.term[j].R[1] = systems[i].solverparams.rho*
            #     #             systems[i].branches.c0[j][end]/systems[i].branches.A0[j];
            #     # end
            # end

            # recalculate and output parameter averages (non-dimensional)
            θhat = mean(θ);
            append!(θout,[θhat])
            # println("Normalized analysis parameter averages: $θhat")

            # analysis parameter variance (for post-processing)
            Ptt = zeros(length(θhat))
            for i = 1:length(soln)
                Ptt += (θ[i] .- θhat).^2;
            end
            Ptt ./= ensemblesize;
            append!(Pθout,[Ptt])

            # # find new time step to satisfy CFL condition for all members
            # h = [];
            # for i = 1:length(soln)
            #     for j in 1:length(systems[i].branches.ID)
            #         push!(h,systems[i].solverparams.CFL*systems[i].branches.k[j]/
            #             systems[i].branches.c0[j][end]);
            #     end
            # end
            # if minimum(h) != systems[1].solverparams.h
            #     # println("Old time step size: $(systems[1].solverparams.h) s")
            #     for i = 1:length(soln)
            #         systems[i].solverparams.h = minimum(h);
            #     end
            #     println("New time step size: $(systems[1].solverparams.h) s")
            #     # re-discretize
            #     if systems[1].solverparams.numbeats < systems[1].solverparams.numbeatstotal
            #         systems = pmap((a1,a2)->LegHemodynamics.rediscretizet!(a1,a2),systems,n);
            #     end
            # end
            #
            # # corrected forecast w/ analysis parameters
            # rflag = "assim";
            # soln = pmap((a1,a2,a3,a4,a5,a6,a7)->LegHemodynamics.advancetime!(a1,a2,a3,a4,a5,a6,a7;
            #     runtype=rflag),systems,spl,times,n,nsamp,term_itr,split_itr);
            # systems = [soln[i][1] for i=1:ensemblesize];
            # n = [soln[i][2] for i=1:ensemblesize];
            # times = [soln[i][3] for i=1:ensemblesize];

            # second forecast state, measurements
            X = [zeros(nstates) for i=1:ensemblesize];
            Y = [zeros(nmeas) for i=1:ensemblesize];

            for i = 1:ensemblesize
                X[i][1] = systems[i].branches.Q[14][1,n[1]]/
                    (systems[i].branches.A[14][1,n[1]]*LegHemodynamics.us);
                X[i][2] = systems[i].branches.Q[24][1,n[1]]/
                    (systems[i].branches.A[24][1,n[1]]*LegHemodynamics.us);
                X[i][3] = systems[i].branches.Q[32][1,n[1]]/
                    (systems[i].branches.A[32][1,n[1]]*LegHemodynamics.us);
                Y[i][1] = systems[i].branches.Q[14][1,n[1]]/
                    (systems[i].branches.A[14][1,n[1]]*LegHemodynamics.us);
                Y[i][2] = systems[i].branches.Q[24][1,n[1]]/
                    (systems[i].branches.A[24][1,n[1]]*LegHemodynamics.us);
                Y[i][3] = systems[i].branches.Q[32][1,n[1]]/
                    (systems[i].branches.A[32][1,n[1]]*LegHemodynamics.us);
            end

            # second forecast mean state, measurement
            xhat = mean(X);
            yhat = mean(Y);
            # println("Normalized ̂x, second forecast: $xhat")
            # println("Normalized ̂y, second forecast: $yhat")

            # output second forecast measurement
            append!(yout,[yhat])
            append!(innov,[yhat-y])

            # output ensemble average state
            append!(xout,[xhat])

            # state variance (for post-processing)
            Pxx = zeros(length(xhat))
            for i = 1:length(soln)
                Pxx += (X[i] .- xhat).^2;
            end
            Pxx ./= ensemblesize;
            append!(Pxout,[Pxx])

            # state 2-sigma quantiles (for post-processing)
            for i = 1:length(xhat)
                xq = Float64[];
                for j = 1:length(soln)
                    push!(xq,X[j][i])
                end
                q = quantile(xq,[0.025,0.975]);
                if i == 1
                    append!(lbx,[ones(1)*q[1]])
                    append!(ubx,[ones(1)*q[2]])
                else
                    push!(lbx[end],q[1])
                    push!(ubx[end],q[2])
                end
            end

            # output measurement time
            push!(tout,systems[1].t[n[1]])
        end
    end
    tt = toc();
    for i = 1:ensemblesize
        times[i].tt += tt;
    end
end

# reshape ensemble averages into vectors of individual time series
if assimflag == "yes"
    for i in 1:length(xout[1]) # i indexes state variables
        xo = Float64[];
        Pxo = Float64[];
        lb = Float64[];
        ub = Float64[];
        for j in 1:length(xout) # j indexes time steps
            push!(xo,xout[j][i])
            push!(Pxo,Pxout[j][i])
            push!(lb,lbx[j][i])
            push!(ub,ubx[j][i])
        end
        append!(xoutv,[xo])
        append!(Pxoutv,[Pxo])
        append!(lbxv,[lb])
        append!(ubxv,[ub])
    end

    for i in 1:length(θout[1]) # i indexes state variables
        θo = Float64[];
        Pθo = Float64[];
        for j in 1:length(θout) # j indexes time steps
            push!(θo,θout[j][i])
            push!(Pθo,Pθout[j][i])
        end
        append!(θoutv,[θo])
        append!(Pθoutv,[Pθo])
    end

    for i in 1:length(yout[1])
        yo = Float64[];
        for j in 1:length(yout)
            push!(yo,yout[j][i])
        end
        append!(youtv,[yo])
    end
end

if saveflag == "yes"
    if assimflag == "no"
        file = MAT.matopen(savefile, "w")
        write(file,"system",system)
        close(file)
    elseif assimflag == "yes"
        for i=1:ensemblesize
            file = MAT.matopen(savefiles[i],"w");
            write(file,"system",systems[i])
            close(file)
        end
        vnames = ["t" "x" "Px" "theta" "thetascale" "Pt" "lb" "ub"];
        for i in 1:length(vnames)
            file = MAT.matopen("$(vnames[i])_converged_small_dev_1.mat","w");
            if vnames[i] == "t"
                write(file,"t",tout)
            elseif vnames[i] == "x"
                write(file,"x",xoutv)
            elseif vnames[i] == "Px"
                write(file,"Px",Pxoutv)
            elseif vnames[i] == "theta"
                write(file,"theta",θoutv)
            elseif vnames[i] == "thetascale"
                write(file,"thetascale",θs)
            elseif vnames[i] == "Pt"
                write(file,"Pt",Pθoutv)
            elseif vnames[i] == "lb"
                write(file,"lb",lbxv)
            elseif vnames[i] == "ub"
                write(file,"ub",ubxv)
            end
            close(file)
        end
    end
end

if assimflag == "no"
    if timeflag == "yes"
        return system, n, times
    elseif timeflag == "no"
        return system, n
    end
elseif assimflag == "yes"
    return n,tout,xoutv,youtv,ioutv,θoutv,Pθoutv,Pxoutv,lbxv,ubxv,θs
end

end
