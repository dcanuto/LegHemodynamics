type ArterialTerminal # lumped peripheral vasculature
    C::Vector{Float64}
    R::Vector{Float64}
    P::Any
    Q::Any

    function ArterialTerminal()
        this = new()
        this.C = Vector{Float64}[];
        this.R = Vector{Float64}[];
        this.P = Array{Float64,2}[];
        this.Q = Array{Float64,2}[];
        return this
    end
end

type ArterialBranches # 1D arterial domain
    name::Vector{String}
    parentname::Vector{String}
    ID::Vector{Int64}
    parentID::Vector{Int64}
    lengthinmm::Vector{Float64}
    radiusinmm::Vector{Float64}
    A0::Any
    beta::Any
    c0::Any
    k::Vector{Float64}
    A::Any
    Q::Any
    P::Any
    Fp::Any
    Fbarforward::Any
    Fbarbackward::Any
    Abackward::Any
    Aforward::Any
    Qbackward::Any
    Qforward::Any
    W1end::Vector{Float64}
    W1::Vector{Float64}
    W2::Vector{Float64}
    W1root::Float64
    W2root::Float64
    children::Any
    term::Any

    function ArterialBranches(filename="test.csv",old=Dict("a"=>0),restart="no")
        this = new()
        this.name = Vector{String}[];
        this.parentname = Vector{String}[];
        this.ID = Vector{Int64}[];
        this.parentID = Vector{Int64}[];
        this.children = Vector{Int64}[];
        this.lengthinmm = Vector{Float64}[];
        this.radiusinmm = Vector{Float64}[];
        # this.A0 = Array{Float64,1}[];
        this.A0 = Vector{Float64}[];
        this.beta = Vector{Float64}[];
        this.c0 = Array{Float64,1}[];
        this.k = Vector{Float64}[];
        this.A = Array{Float64,2}[];
        this.Q = Array{Float64,2}[];
        this.P = Array{Float64,2}[];
        this.Fp = Array{Float64,1}[];
        this.Fbarforward = Array{Float64,1}[];
        this.Fbarbackward = Array{Float64,1}[];
        this.Abackward = Array{Float64,1}[];
        this.Aforward = Array{Float64,1}[];
        this.Qbackward = Array{Float64,1}[];
        this.Qforward = Array{Float64,1}[];
        this.W1end = Vector{Float64}[];
        this.W1 = Vector{Float64}[];
        this.W2 = Vector{Float64}[];
        if restart == "no"
            this.W1root = Float64(0);
            this.W2root = Float64(0);

            # pull in artery data from text file
            temp = LegHemodynamics.loadtexttree(filename);
            this.name = [string(temp[1,:Name])]
            this.parentname = [string(temp[1,:ParentName])]
            this.ID = [temp[1,:ID]]
            this.parentID = [temp[1,:parentID]]
            this.lengthinmm = [temp[1,:length_mm]]
            this.radiusinmm = [temp[1,:radius_mm]]
            factor = 1.0;
            this.A0 = [factor*temp[1,:A0_m2]]
            this.beta = [temp[1,:beta_Pam]/factor]
            for i = 2:length(temp[:Name])
                push!(this.name,string(temp[i,:Name]))
                push!(this.parentname,string(temp[i,:ParentName]))
                push!(this.ID,temp[i,:ID])
                push!(this.parentID,temp[i,:parentID])
                push!(this.lengthinmm,temp[i,:length_mm])
                push!(this.radiusinmm,temp[i,:radius_mm])
                push!(this.A0,factor*temp[i,:A0_m2])
                push!(this.beta,temp[i,:beta_Pam]/factor)
            end
            for i in 1:size(temp[:children_1],1)
                if !isa(temp[i,:children_1],Missings.Missing)
                    push!(this.children,[temp[i,:children_1]])
                else
                    push!(this.children,[])
                end
                if !isa(temp[i,:children_2],Missings.Missing)
                    push!(this.children[i],temp[i,:children_2])
                end
                if !isa(temp[i,:children_3],Missings.Missing)
                    push!(this.children[i],temp[i,:children_3])
                end
            end
        elseif restart == "yes"
            temp = LegHemodynamics.loadtexttree("arterylist_2.txt");
            this.name = [string(temp[1,:Name])]
            this.parentname = [string(temp[1,:ParentName])]
            this.W1root = old["W1root"];
            this.W2root = old["W2root"];
            # this.name = [old["name"][1]];
            # this.parentname = [old["parentname"][1]];
            this.ID = [old["ID"][1]];
            this.parentID = [old["parentID"][1]];
            this.lengthinmm = [old["lengthinmm"][1]];
            this.radiusinmm = [old["radiusinmm"][1]];
            # this.A0 = [old["A0"][1]];
            this.beta = [old["beta"][1]]
            for i = 2:length(old["ID"])
                push!(this.name,string(temp[i,:Name]))
                push!(this.parentname,string(temp[i,:ParentName]))
                # push!(this.name,old["name"][i]);
                # push!(this.parentname,old["parentname"][i]);
                push!(this.ID,old["ID"][i]);
                push!(this.parentID,old["parentID"][i]);
                push!(this.lengthinmm,old["lengthinmm"][i]);
                push!(this.radiusinmm,old["radiusinmm"][i]);
                # push!(this.A0,old["A0"][i]);
                push!(this.beta,old["beta"][i]);
            end
            for i = 1:length(old["children"])
                temp = Array{Int64}(1);
                if length(old["children"][i]) == 1
                    temp = [old["children"][i]];
                elseif length(old["children"][i]) > 1
                    temp = old["children"][i];
                else
                    temp = [0];
                end
                push!(this.children,temp)
                deleteat!(this.children[i],findin(this.children[i],0))
            end
        end

        this.term = Vector{LegHemodynamics.ArterialTerminal}(length(this.ID));

        return this
    end
end
