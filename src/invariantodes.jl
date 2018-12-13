function invariantodes!(dW::Vector{Float64},W1::Vector{Float64},W2::Vector{Float64},
    lf::Float64,lb::Float64,Q::Vector{Float64},A::Vector{Float64},
    rho::Float64,μ::Float64,diff::Float64,k::Float64,ts::Float64,xs::Float64,
    Qs::Float64,As::Float64)

    dW[1] = -ts.*Qs./(xs.*As).*lf.*dot([-1. 1.]./k,W1[1:2]) .-
        ts./As.*diff.*π.*μ./rho.*Q[end]./(A[end].^2);
    dW[2] = -ts.*Qs./(xs.*As).*lb.*dot([-1. 1.]./k,W2[1:2]) .-
        ts./As.*diff.*π.*μ./rho.*Q[1]./(A[1].^2);

end
