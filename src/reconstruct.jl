function reconstruct(w::Vector{Float64},c::Vector{Float64},U::Vector{Float64},
    j::Int,k::Float64,σ::Vector{Float64})
    P = function(Pout::Vector{Float64},x::Float64,Pi::Vector{Float64})
        if j == 1 # left boundary polynomial
            σ[1] = U[2] .- U[1];
            σ[2] = U[3] .- U[2];
            Pi[1] = U[1];
            Pi[2] = U[2] .+ σ[1]./k.*(x .- 1.5.*k);
            Pi[3] = 1./c[3].*(-c[1].*U[1] .+ (1-c[2]).*U[2] .- 1./24.*(σ[2] - σ[1]) .+
                (σ[2] .+ σ[1] .- 2.*c[2].*σ[1])./(2.*k).*(x .- 1.5.*k) .+ (σ[2] - σ[1])./(2.*k.^2).*
                (x .- 1.5.*k).^2);
        elseif j == length(U) # right boundary polynomial
            σ[1] = U[end-2] .- U[end-1];
            σ[2] = U[end-1] .- U[end];
            Pi[1] = U[end];
            Pi[2] = U[end-1] .- σ[2]./k.*(x - k.*(length(U)-1.5));
            Pi[3] = 1./c[3].*(-c[1].*U[end] .+ (1-c[2]).*U[end-1] .- 1./24.*(σ[1] - σ[2]) .-
                ((σ[1] .+ σ[2]) - 2.*c[2].*σ[2])./(2.*k).*(x .- k.*(length(U)-1.5)) .+
                (σ[1] - σ[2])./(2.*k.^2).*(x .- k.*(length(U)-1.5)).^2);
        else # interior polynomial
            σ[1] = U[j] .- U[j-1];
            σ[2] = U[j+1] .- U[j];
            Pi[1] = U[j] .+ σ[1]./k.*(x .- k.*(j-0.5));
            Pi[2] = 1./c[2].*((1.-c[1].-c[3]).*U[j] .- 1./24.*(σ[2] .- σ[1]) .+
                ((1 - 2.*c[3]).*U[j+1] .+ 2.*(c[3] .- c[1]).*U[j] .+
                (2.*c[1] .- 1).*U[j-1])./(2.*k).*(x .- k.*(j-0.5)) .+
                (σ[2] .- σ[1])./(2.*k.^2).*(x .- k.*(j-0.5)).^2);
            Pi[3] = U[j] .+ σ[2]./k.*(x .- k.*(j-0.5));
        end
        Pout[1] = dot(w,Pi);
    end
    return P
end
