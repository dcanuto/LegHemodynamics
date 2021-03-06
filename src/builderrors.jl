type Errors
    a::Float64 # parameter distribution smoothing
    h::Float64

    odev::Vector{Float64} # σ for observation, param distributions
    pdev::Vector{Float64}

    lb::Vector{Float64} # lower/upper bounds for truncated normal param distrs.
    ub::Vector{Float64}

    # Amin::Vector{Float64} # minimum values for A0 for non-measured arteries
    # bmax::Vector{Float64}

    function Errors(nparams=11)
        this = new();
        δ = 0.99; # ∃ (0, 1, typically within 0.95-0.99), lower values = less parameter dispersion
        this.a = (3*δ-1)/(2*δ);
        this.h = sqrt.(1-this.a^2);
        this.odev = [25*mmTom,25*mmTom,25*mmTom];
        # this.pdev = 1e-12*ones(nparams);
        this.pdev = [1e9,2e9,1e10,4e-11,2e-11,2e-12,0.2,0.2,0.2,0.2,0.2,0.2,8e6,8e6,1e7];
        # this.pdev = [1e9,2.3e9,1e10,4e-11,2e-11,4.8e-12,0.1,0.3,0.05,0.05,0.05,4e-7,4.3e-7,2e-7];
        # this.pdev = [1e9,2.3e9,1e10,4e-11,2e-11,4.8e-12,0.1,0.3,7e5,7e5,7e5];
        # this.pdev = [1e9,2.3e9,1e10,4e-11,2e-11,4.8e-12,0.1,0.1,4e-7,4.3e-7,2e-7];
        this.lb = -Inf*ones(nparams);
        this.ub = Inf*ones(nparams);
        this.lb[1] = 1e10; # R1
        this.lb[2] = 1e10; # R2
        this.lb[3] = 1e10; # R3
        this.lb[4] = 1e-10; # C1
        this.lb[5] = 1e-10; # C2
        this.lb[6] = 7e-12; # C3
        this.lb[7] = 0.1; # Cfac
        this.lb[8] = 0.1; # Rfac
        this.lb[9] = 0.1; # Cfac
        this.lb[10] = 0.1; # Rfac
        this.lb[11] = 0.1; # Cfac
        this.lb[12] = 0.1; # Rfac
        # this.lb[9] = 2e-6; # A1
        # this.lb[10] = 2.15e-6; # A2
        # this.lb[11] = 9.5e-7; # A3
        this.lb[13] = 2e7; # β1
        this.lb[14] = 2e7; # β2
        this.lb[15] = 2e7; # β3
        this.ub[1] = 1e11;
        this.ub[2] = 1e11;
        this.ub[3] = 1e11;
        this.ub[4] = 1e-9;
        this.ub[5] = 1e-9; #2.5e-10
        this.ub[6] = 2e-11;
        # this.ub[7] = 1.5;
        # this.ub[8] = 1.2;
        # this.ub[9] = 1.5;
        # this.ub[10] = 1.2;
        # this.ub[11] = 1.5;
        # this.ub[12] = 1.2;
        this.ub[13] = 2e8; # β1 8.5e7
        this.ub[14] = 2e8; # β2 1.5e8
        this.ub[15] = 2e8; # β3 5.5e8
        # this.ub[9] = 6e7; # E1
        # this.ub[10] = 6e7; # E2
        # this.ub[11] = 6e7; # E3

        # pull in artery data from text file
        temp = LegHemodynamics.loadtexttree("arterylist_2.txt");
        # this.Amin = 0.8*[temp[1,:A0_m2]];
        # for i = 2:length(temp[:A0_m2])
        #     push!(this.Amin,0.8*temp[i,:A0_m2])
        # end
        # this.bmax = 3*[temp[1,:beta_Pam]];
        # for i = 2:length(temp[:beta_Pam])
        #     push!(this.bmax,3*temp[i,:beta_Pam])
        # end
        return this
    end
end
