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
    # terminal properties from Joe's data
    # Rdefault = [74.7761,136.1894,110.4835,98.8103,170.0850,247.8310,16.8961,
    #     103.3964,280.8708,184.7235,143.7461,34.2842,51.0196,144.4766,250.2162,
    #     138.8846,92.7085]*gTokg/(mmTom^4);
    Rdefault = 0.75*[77.6909,141.4980,114.7901,102.6618,176.7149,257.4914,18.6383,
        107.4268,291.8191,191.9240,149.3493,35.9008,53.0083,150.1083,259.9696,
        142.5643,95.1046]*gTokg/(mmTom^4);
    Cdefault = [0.08972235,0.04926296,0.06072485,0.06789874,0.03944552,0.02707124,
        0.39707822,0.06488707,0.02388675,0.03631963,0.04667319,0.19569034,0.13150028,
        0.04643721,0.02681318,0.04830696,0.00723676]*(mmTom^4)/gTokg;
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
