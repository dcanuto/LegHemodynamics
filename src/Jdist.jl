function Jdist!(J::Vector{Float64},x::Vector{Float64},C::Float64,R0::Float64,
    R1::Float64,Ws::Float64,W1::Float64,rho::Float64,beta::Float64,c0::Float64,
    a1::Float64,a2::Float64,h::Float64)
    J[1] = -a1/(4*a2*Ws^4)*(0.125*(W1-Ws*x[1])+c0) - (-3*W1^4/(8192*Ws^4) +
        W1^3*x[1]^2/(2048*Ws^3) - W1^3*c0/(128*Ws^4) + 3*W1^2*x[1]^2/(4096*Ws^2) -
        3*W1^2*c0^2/(64*Ws^4) - 3*W1*x[1]^3/(2048*Ws) + 3*W1*x[1]^2*c0/(128*Ws^2) -
        3*W1*x[1]*c0^2/(32*Ws^3) + 5*x[1]^4/8192 - x[1]^3*c0/(64*Ws) +
        9*x[1]^2*c0^2/(64*Ws^2) - x[1]*c0^3/(2*Ws^3) + c0^4/(2*Ws^4));
end
