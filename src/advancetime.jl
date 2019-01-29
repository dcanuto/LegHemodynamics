function advancetime!(system::LegSystem,spl::Dierckx.Spline1D,times::CVTimer,
    n::Int64,nsamp::Int64,terms::Vector{Int64},splits::Vector{Int64};runtype="predict")
    tic();
    if runtype == "predict"
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
            times.tr += toq();
            if mod(n/nsamp,1) == 0
                n+=1;
                break
            end
            n+=1;
        end
    elseif runtype == "assim"
        # check if new cardiac cycle has just begun
        numbeats = system.solverparams.numbeats;
        ttotal = system.solverparams.th*(system.solverparams.numbeats+1);
        if system.t[n+1] > ttotal
            println("New cardiac cycle beginning. Decrementing number of cycles for assimilation.")
            numbeats -= 1;
        end
        # time within heart cycle
        tp = system.t[n] - sum(system.solverparams.th*numbeats);
        # interpolate proximal flow velocity from patient data
        uprox = spl(tp);
        # TVD RK3 time integration
        LegHemodynamics.tvdrk3!(system,times,n-1,splits,terms,-uprox);
        tic();
        # arterial pressure update
        LegHemodynamics.arterialpressure!(system,n-1);
        times.tr += toq();
    end
    if runtype == "predict"
        times.tt += toc();
    elseif runtype == "assim"
        times.tt = toq();
    end
    return system,n,times
end
