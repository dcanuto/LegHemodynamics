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
    system = LegHemodynamics.buildall(loadfile;numbeatstotal=15,restart=rstflag);
    savefile = "test.mat" # filename for saving (only used if saveflag == "yes")
elseif assimflag == "yes"
    ensemblesize = 3;
    if rstflag == "no"
        loadfiles = ["arterylist_2.txt" for i=1:ensemblesize];
    elseif rstflag == "yes"
        loadfiles = ["test_1_$i.mat" for i=1:ensemblesize];
    end
    systems = pmap((a1)->LegHemodynamics.buildall(a1;numbeatstotal=3,restart=rstflag),loadfiles);
    savefiles = ["test_1_$i.mat" for i=1:ensemblesize];
end

# create interpolation object for upstream flow velocity data (needed for proximal BC)
patientdata = MAT.matread("flowrates_mL.mat");
tdata = reshape(patientdata["t"],101);
udata = reshape(patientdata["u_pop_cms"],101)*LegHemodynamics.cmTom;
if assimflag == "no"
    spl = Dierckx.Spline1D(tdata,udata;k=1);
elseif assimflag == "yes"
    spl = [Dierckx.Spline1D(tdata,udata;k=1) for i=1:ensemblesize];
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
    nsamp = [1000 for i=1:ensemblesize];
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
            println("Current time step: $(n[1]-1)")
            println("Current time: $(systems[1].t[n[1]])")
        end
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
    end
end

if assimflag == "no"
    if timeflag == "yes"
        return system, n, times
    elseif timeflag == "no"
        return system, n
    end
elseif assimflag == "yes"
    return n
end

end
