function abseigs1d(abso::Vector{Float64},eigs::Vector{Float64},J::Matrix{Float64},
    Q::Float64,A::Float64,beta::Float64,rho::Float64,ts::Float64,xs::Float64,
    Qs::Float64,As::Float64)
    J1d(J,Q,A,beta,rho,ts,xs,Qs,As);
    eigs[1] = 0.5*(J[1,1] + J[2,2]) + sqrt.(0.25*(J[1,1] + J[2,2])^2 -
        J[1,1]*J[2,2] + J[1,2]*J[2,1]);
    eigs[2] = 0.5*(J[1,1] + J[2,2]) - sqrt.(0.25*(J[1,1] + J[2,2])^2 -
        J[1,1]*J[2,2] + J[1,2]*J[2,1]);
    abso[1] = abs(eigs[1]);
    abso[2] = abs(eigs[2]);
end
