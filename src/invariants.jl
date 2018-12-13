function invariants!(W1::Vector{Float64},W2::Vector{Float64},lf::Vector{Float64},
    lb::Vector{Float64},A::Vector{Float64},Q::Vector{Float64},
    c0::Float64,beta::Float64,rho::Float64,μ::Float64,diff::Float64,
    Qs::Float64,As::Float64)

    # characteristic speeds
    lb[1] = (Q[1]./A[1].*Qs./As .- sqrt.(0.5.*beta./rho).*(A[1].*As).^0.25).*As./Qs;
    lf[1] = (Q[end]./A[end].*Qs./As .+ sqrt.(0.5.*beta./rho).*(A[end].*As).^0.25).*As./Qs;

    # invariants
    W1[:] .= (Q[end-1:end]./A[end-1:end].*Qs./As .+
        4.*(sqrt.(0.5.*beta./rho).*(A[end-1:end].*As).^0.25 .- c0)).*As./Qs;
    W2[:] .= (Q[1:2]./A[1:2].*Qs./As .-
        4.*(sqrt.(0.5.*beta./rho).*(A[1:2].*As).^0.25 .- c0)).*As./Qs;

end
