type LegSystem # entire solution
    branches::LegHemodynamics.ArterialBranches
    solverparams::LegHemodynamics.SolverParams
    t::Vector{Float64}

    function LegSystem(filename="test.csv",restart="no")
        this = new()
        if restart == "no"
            this.branches = LegHemodynamics.ArterialBranches(filename);
        elseif restart == "yes"
            vars = MAT.matread(filename);
            sys = vars["system"];
            branches = sys["branches"];
            this.branches = LegHemodynamics.ArterialBranches(filename,branches,restart)
        else
            error("Keyword restart must either be yes or no. Aborting.")
        end
        this.solverparams = LegHemodynamics.SolverParams();
        this.t = Vector{Float64}[];
        return this
    end
end

# build solution struct
function buildall(filename="test.csv";numbeatstotal=1,restart="no",injury="no")
    # # terminal properties adapted from Danielsen (1998)
    # Rdefault = [0.3,0.21,0.003,0.01]*mmHgToPa/cm3Tom3;
    # Cdefault = [0.01,1.64,1.81,13.24,73.88]*cm3Tom3/mmHgToPa;
    # Vdefault = [370.,370.,400.,500.,1400.]*cm3Tom3;
    # Ldefault = 5e-5*mmHgToPa/cm3Tom3;
    if restart == "no"
        system = LegHemodynamics.LegSystem(filename);
        system.solverparams.numbeatstotal = numbeatstotal;
        LegHemodynamics.calcbranchprops!(system);
        LegHemodynamics.discretizebranches!(system);
        # LegHemodynamics.assignterminals!(system,Rdefault,Cdefault,Vdefault,Ldefault,lowerflowfraction,
        #     venousfractionofR,venousfractionofL,venousfractionofC,venousfractionofV0);
        # LegHemodynamics.discretizeperiphery!(system);
        # LegHemodynamics.applybranchics!(system);
        # LegHemodynamics.applyperipheryics!(system);
        # LegHemodynamics.applycustomics!(system);
    elseif restart == "yes"
        vars = MAT.matread(filename);
        sys = vars["system"];
        branches = sys["branches"];
        term = branches["term"];
        system = LegHemodynamics.LegSystem(filename,restart);
        system.solverparams.numbeatstotal = numbeatstotal;
        LegHemodynamics.calcbranchprops!(system,branches,restart);
        LegHemodynamics.discretizebranches!(system,sys,restart);
        LegHemodynamics.assignterminals!(system,Rdefault,Cdefault,Vdefault,Ldefault,lowerflowfraction,
            venousfractionofR,venousfractionofL,venousfractionofC,venousfractionofV0,term,restart);
        LegHemodynamics.discretizeperiphery!(system);
        LegHemodynamics.applybranchics!(system,sys,restart);
        LegHemodynamics.applyperipheryics!(system,sys,restart);
    end
    # LegHemodynamics.updatevolumes!(system,0);
    return system
end
