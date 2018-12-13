function weights!(w::Vector{Float64},c::Vector{Float64},ϵ::Float64,IS::Vector{Float64},
    s::Vector{Float64})
    for i = 1:length(w)
        w[i] = c[i]/((ϵ + IS[i])^2);
    end
    s[1] = sum(w);
    for i = 1:length(w)
        w[i] /= s[1];
    end
end
