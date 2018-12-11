function discretizebranches!(system::CVSystem,old=Dict("a"=>0),restart="no")
    h = [];

    # branch grid spacing
    for i in 1:length(system.branches.ID)
        push!(system.branches.k,system.branches.lengthincm[i]*cmTom/
            (system.solverparams.JL-1));
        push!(h,system.solverparams.CFL*system.branches.k[i]/
            system.branches.c0[i][end]);
    end

    # time step guaranteed to satisfy CFL for all branches
    # system.solverparams.h = minimum(h);

    # fix time step size to couple with 3D liver model
    system.solverparams.h = 1e-4;

    if restart == "no"
        system.solverparams.numsteps = ceil(system.heart.activation.th[1]/
            system.solverparams.h);

        # allocate space for 1D domain solution variables
        for i in 1:length(system.branches.ID)
            push!(system.branches.A,
                zeros(system.solverparams.JL,system.solverparams.numsteps+1));
            push!(system.branches.Q,
                zeros(system.solverparams.JL,system.solverparams.numsteps+1));
            push!(system.branches.P,
                zeros(system.solverparams.JL,system.solverparams.numsteps+1));
            push!(system.branches.Fp,
                zeros(2*system.solverparams.JL));
            push!(system.branches.Fbarforward,
                zeros(2*system.solverparams.JL-4));
            push!(system.branches.Fbarbackward,
                zeros(2*system.solverparams.JL-4));
            push!(system.branches.Abackward,
                zeros(system.solverparams.JL-2));
            push!(system.branches.Aforward,
                zeros(system.solverparams.JL-2));
            push!(system.branches.Qbackward,
                zeros(system.solverparams.JL-2));
            push!(system.branches.Qforward,
                zeros(system.solverparams.JL-2));
            push!(system.branches.W1end,0.);
            push!(system.branches.W1,0.);
            push!(system.branches.W2,0.);
        end

        # discretize time for first cardiac cycle
        system.t = system.solverparams.h*[0:1:size(system.branches.A[1],2)-1;];
    elseif restart == "yes"
        # determine time shift
        temp = old["heart"];
        temp = temp["activation"];
        system.solverparams.tshift = old["t"][end] - sum(temp["th"][1:end-1]);

        # change grid spacing (e.g., for grid refinement study)
        JLnew = 11;

        if JLnew != system.solverparams.JL
            system.solverparams.JL = JLnew;
            for i = 1:length(system.branches.ID)
                system.branches.k[i] = system.branches.lengthincm[i]*cmTom/
                    (system.solverparams.JL-1);
                h[i] = system.solverparams.CFL*system.branches.k[i]/
                    system.branches.c0[i][end];
            end

            # fix other index ranges
            system.solverparams.acols = [1:system.solverparams.JL;];
            system.solverparams.qcols = [system.solverparams.JL+1:2*system.solverparams.JL;];
            system.solverparams.acolspre = [2:system.solverparams.JL-1;];
            system.solverparams.qcolspre = [system.solverparams.JL+2:2*system.solverparams.JL-1;];
            system.solverparams.acolscor = [1:system.solverparams.JL-2;];
            system.solverparams.qcolscor = [system.solverparams.JL-1:2*system.solverparams.JL-4;];
            system.solverparams.colsint = [2:system.solverparams.JL-1;];

            # time step guaranteed to satisfy CFL for all branches
            system.solverparams.h = minimum(h);
        end

        # update total number of time steps
        ntoadd = ceil((system.heart.activation.th[end]-
            system.solverparams.tshift)/system.solverparams.h);
        system.solverparams.numsteps=ntoadd;

        # update discrete times
        ttoadd = [0:1:ntoadd;]*system.solverparams.h + system.solverparams.tshift;
        system.t = ttoadd;

        # allocate space for 1D domain solution variables
        for i in 1:length(system.branches.ID)
            push!(system.branches.A,
                zeros(system.solverparams.JL,system.solverparams.numsteps+1));
            push!(system.branches.Q,
                zeros(system.solverparams.JL,system.solverparams.numsteps+1));
            push!(system.branches.P,
                zeros(system.solverparams.JL,system.solverparams.numsteps+1));
            push!(system.branches.Fp,
                zeros(2*system.solverparams.JL));
            push!(system.branches.Fbarforward,
                zeros(2*system.solverparams.JL-4));
            push!(system.branches.Fbarbackward,
                zeros(2*system.solverparams.JL-4));
            push!(system.branches.Abackward,
                zeros(system.solverparams.JL-2));
            push!(system.branches.Aforward,
                zeros(system.solverparams.JL-2));
            push!(system.branches.Qbackward,
                zeros(system.solverparams.JL-2));
            push!(system.branches.Qforward,
                zeros(system.solverparams.JL-2));
            push!(system.branches.W1end,0.);
            push!(system.branches.W1,0.);
            push!(system.branches.W2,0.);
        end
    end
    print("Number of nodes per branch: ")
    println(system.solverparams.JL)
    return system
end
