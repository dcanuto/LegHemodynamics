type Errors
    a::Float64 # parameter distribution smoothing
    h::Float64

    odev::Vector{Float64} # σ for observation, param distributions
    pdev::Vector{Float64}

    lb::Vector{Float64} # lower/upper bounds for truncated normal param distrs.
    ub::Vector{Float64}

    function Errors(nparams=3)
        this = new();
        δ = 0.99; # ∃ (0, 1, typically within 0.95-0.99), lower values = less parameter dispersion
        this.a = (3*δ-1)/(2*δ);
        this.h = sqrt.(1-this.a^2);
        this.odev = [20*mmTom,20*mmTom,20*mmTom];
        # this.pdev = 1e-12*ones(nparams);
        this.pdev = [1e9,2.3e9,1e10,4e-11,2e-11,4.8e-12];
        this.lb = -Inf*ones(nparams);
        this.ub = Inf*ones(nparams);
        this.lb[1] = 1e9; # R1
        this.lb[2] = 2.3e9; # R2
        this.lb[3] = 1e10; # R3
        this.lb[4] = 8e-11; # C1
        this.lb[5] = 4e-11; # C2
        this.lb[6] = 1e-11; # C3
        # this.lb[7] = 1.; # m1
        # this.lb[8] = 1; # m2
        # this.lb[9] = 1e6; # Emax
        # this.ub[5] = 0.4; # τ1
        # this.ub[6] = 0.55; # Δτ
        # this.ub[9] = 2e8; # Emax
        return this
    end
end
