function fdist!(f::Vector{Float64},x::Vector{Float64},D::Vector{Float64},C::Float64,
    R0::Float64,R1::Float64,Ws::Float64,W1::Float64,rho::Float64,beta::Float64,
    c0::Float64,a1::Float64,a2::Float64,Pn::Float64,h::Float64)
    f[1] = D[1]*(a1/(a2*Ws^5)*(0.125*(W1-Ws*x[1])+c0)^2 - (W1^5/(8192*Ws^5) -
        3*W1^4*x[1]/(8192*Ws^4) + W1^4*c0/(256*Ws^5) + W1^3*x[1]^2/(4096*Ws^3) -
        W1^3*x[1]*c0/(128*Ws^4) + 3*W1^3*c0^2/(64*Ws^5) + W1^2*x[1]^3/(4096*Ws^2) -
        3*W1^2*x[1]*c0^2/(64*Ws^4) + W1^2*c0^3/(4*Ws^5) - 3*W1*x[1]^4/(8192*Ws) +
        W1*x[1]^3*c0/(128*Ws^2) - 3*W1*x[1]^2*c0^2/(64*Ws^3) + W1*c0^4/(2*Ws^5) +
        x[1]^5/8192 - x[1]^4*c0/(256*Ws) + 3*x[1]^3*c0^2/(64*Ws^2) -
        x[1]^2*c0^3/(4*Ws^3) + x[1]*c0^4/(2*Ws^4)) + Pn/(a2*Ws^5));
end
