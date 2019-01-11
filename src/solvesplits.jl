function solvesplits!(iters::Vector{Int64},Q::Vector{Float64},A::Vector{Float64},
    children::Vector{Int64},Qold::Vector{Float64},Aold::Vector{Float64},W::Vector{Float64},
    beta::Vector{Float64},A0::Vector{Float64},c0::Vector{Float64},rho::Float64,
    f::Function,J::Function,maxiter::Int16,maxval::Float64,epsJ::Float64,epsN::Float64)
    LegHemodynamics.newton!(iters,Q,A,children,Qold,Aold,W,beta,A0,c0,rho,f,J,
        maxiter,maxval,epsJ,epsN);
end
