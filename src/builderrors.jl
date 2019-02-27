type Errors
    a::Float64 # parameter distribution smoothing
    h::Float64

    odev::Vector{Float64} # σ for observation, param distributions
    pdev::Vector{Float64}

    lb::Vector{Float64} # lower/upper bounds for truncated normal param distrs.
    ub::Vector{Float64}

    function Errors(nparams=11)
        this = new();
        δ = 0.99; # ∃ (0, 1, typically within 0.95-0.99), lower values = less parameter dispersion
        this.a = (3*δ-1)/(2*δ);
        this.h = sqrt.(1-this.a^2);
        this.odev = [40*mmTom,40*mmTom,40*mmTom];
        # this.pdev = 1e-12*ones(nparams);
        this.pdev = [1e9,2.3e9,1e10,4e-11,2e-11,4.8e-12,0.3,0.3,4e7,3.8e7,3.8e7];
        # this.pdev = [1e9,2.3e9,1e10,4e-11,2e-11,4.8e-12,0.3,1e-12,4e7,3.8e7,3.8e7];
        # this.pdev = [1e9,4e-11];
        this.lb = -Inf*ones(nparams);
        this.ub = Inf*ones(nparams);
        this.lb[1] = 1e9; # R1
        this.lb[2] = 2.3e9; # R2
        this.lb[3] = 1e10; # R3
        this.lb[4] = 8e-11; # C1
        this.lb[5] = 4e-11; # C2
        this.lb[6] = 1e-11; # C3
        this.lb[7] = 0.1; # Cfac
        this.lb[8] = 0.1; # Rfac
        this.lb[9] = 2e7; # β1
        this.lb[10] = 2e7; # β2
        this.lb[11] = 2e7; # β3
        this.ub[8] = 3;
        return this
    end
end
