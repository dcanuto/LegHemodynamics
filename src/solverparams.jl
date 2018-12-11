type SolverParams
    numbeats::Int64
    numbeatstotal::Int64
    numsteps::Int64
    nstart::Int64
    JL::Int64
    acols::Vector{Int8}
    qcols::Vector{Int8}
    acolspre::Vector{Int8}
    qcolspre::Vector{Int8}
    acolscor::Vector{Int8}
    qcolscor::Vector{Int8}
    colsint::Vector{Int8}
    CFL::Float64
    h::Float64
    tshift::Float64
    rho::Float64
    mu::Float64
    nu::Float64
    diffusioncoeff::Float64
    abserror::Float64
    maxiter::Int16
    epsJ::Float64
    epsN::Float64
    maxval::Float64
    totaliter::Int64

    function SolverParams()
        this = new()
        this.CFL = 0.5;
        this.JL = 11;
        this.acols = [1:this.JL;];
        this.qcols = [this.JL+1:2*this.JL;];
        this.acolspre = [2:this.JL-1;];
        this.qcolspre = [this.JL+2:2*this.JL-1;];
        this.acolscor = [1:this.JL-2;];
        this.qcolscor = [this.JL-1:2*this.JL-4;];
        this.colsint = [2:this.JL-1;];
        this.nstart = 0;
        this.numbeats = 0;
        this.numbeatstotal = 0;
        this.tshift = 0;
        this.rho = 1060;
        this.mu = 0.004;
        this.nu = 0.5;
        this.diffusioncoeff = 22.0;
        this.abserror = 1e-10;
        this.maxiter = 200;
        this.epsJ = 1e-11;
        this.epsN = 1e-8;
        this.maxval = 1e7;
        this.totaliter = 0;
        return this
    end
end
