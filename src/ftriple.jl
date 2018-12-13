function ftriple!(f::Vector{Float64},x::Vector{Float64},beta::Vector{Float64},A0::Vector{Float64},
    rho::Float64,c0::Vector{Float64},W::Vector{Float64})
    # f = zeros(8);
    # conservation of mass
    f[1] = x[1] .- x[3] .- x[5] .- x[7];

    # total P
    f[2] = (beta[1].*(sqrt.(x[2]) .- sqrt.(A0[1])) + 0.5.*rho.*(x[1]/x[2]).^2) .-
        (beta[2].*(sqrt.(x[4]) .- sqrt.(A0[2])) + 0.5.*rho.*(x[3]/x[4]).^2);
    f[3] = (beta[1].*(sqrt.(x[2]) .- sqrt.(A0[1])) + 0.5.*rho.*(x[1]/x[2]).^2) .-
        (beta[3].*(sqrt.(x[6]) .- sqrt.(A0[3])) + 0.5.*rho.*(x[5]/x[6]).^2);
    f[4] = (beta[1].*(sqrt.(x[2]) .- sqrt.(A0[1])) + 0.5.*rho.*(x[1]/x[2]).^2) .-
        (beta[4].*(sqrt.(x[8]) .- sqrt.(A0[4])) + 0.5.*rho.*(x[7]/x[8]).^2);

    # Riemann invariants
    f[5] = (x[1]./x[2]) .+ 4.*(sqrt.(0.5.*beta[1]./rho).*x[2].^0.25 .- c0[1]) .- W[1];
    f[6] = (x[3]./x[4]) .- 4.*(sqrt.(0.5.*beta[2]./rho).*x[4].^0.25 .- c0[2]) .- W[2];
    f[7] = (x[5]./x[6]) .- 4.*(sqrt.(0.5.*beta[3]./rho).*x[6].^0.25 .- c0[3]) .- W[3];
    f[8] = (x[7]./x[8]) .- 4.*(sqrt.(0.5.*beta[4]./rho).*x[8].^0.25 .- c0[4]) .- W[4];

    # return f
end
