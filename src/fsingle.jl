function fsingle!(f::Vector{Float64},x::Vector{Float64},beta::Vector{Float64},A0::Vector{Float64},
    rho::Float64,c0::Vector{Float64},W::Vector{Float64})
    # f = zeros(4);
    # conservation of mass
    f[1] = x[1] .- x[3];

    # total P
    f[2] = (beta[1].*(sqrt.(x[2]) .- sqrt.(A0[1])) + 0.5.*rho.*(x[1]/x[2]).^2) .-
        (beta[2].*(sqrt.(x[4]) .- sqrt.(A0[2])) + 0.5.*rho.*(x[3]/x[4]).^2);

    # Riemann invariants
    f[3] = (x[1]./x[2]) .+ 4.*(sqrt.(0.5.*beta[1]./rho).*x[2].^0.25 .- c0[1]) .- W[1];
    f[4] = (x[3]./x[4]) .- 4.*(sqrt.(0.5.*beta[2]./rho).*x[4].^0.25 .- c0[2]) .- W[2];

    # return f
end
