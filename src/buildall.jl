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
    Rdefault = ones(17)*0.21*mmHgToPa/cm3Tom3;
    Cdefault = ones(17)*0.01*cm3Tom3/mmHgToPa;
    if restart == "no"
        system = LegHemodynamics.LegSystem(filename);
        system.solverparams.numbeatstotal = numbeatstotal;
        LegHemodynamics.calcbranchprops!(system);
        LegHemodynamics.discretizebranches!(system);
        LegHemodynamics.assignterminals!(system,Rdefault,Cdefault);
        LegHemodynamics.discretizeperiphery!(system);
        LegHemodynamics.applybranchics!(system);
        LegHemodynamics.applyperipheryics!(system);
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
        LegHemodynamics.assignterminals!(system,Rdefault,Cdefault,term,restart);
        LegHemodynamics.discretizeperiphery!(system);
        LegHemodynamics.applybranchics!(system,sys,restart);
        LegHemodynamics.applyperipheryics!(system,sys,restart);
    end
    return system
end
