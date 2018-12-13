function J1d(J::Matrix{Float64},Q::Float64,A::Float64,beta::Float64,rho::Float64,ts::Float64,
    xs::Float64,Qs::Float64,As::Float64)
    if A < 0.
        error("A negative. A = $A.")
    end
    J[1,2] = ts*Qs/(As*xs);
    J[2,1] = -Qs*ts/(As*xs)*Q^2/(A^2) + As^1.5*ts/(Qs*xs)*beta*sqrt(A)/(2.*rho);
    J[2,2] = Qs*ts/(As*xs)*2.*Q/A;
end
