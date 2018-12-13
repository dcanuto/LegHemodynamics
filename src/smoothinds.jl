function smoothinds!(IS::Vector{Float64},U::Vector{Float64},c::Vector{Float64},
    k::Float64,j::Int,σ::Vector{Float64})
    if j == 1 # left boundary
        σ[1] = U[2] - U[1];
        σ[2] = U[3] - U[2];
        IS[2] = σ[1]^2;
        IS[3] = 1/(c[3]^2)*(4/3*σ[2]^2 + (c[2] - 11/3)*σ[2]*σ[1] +
            (10/3 - 3*c[2] + c[2]^2)*σ[1]^2);
    elseif j == length(U) # right boundary
        σ[1] = U[end-2] - U[end-1];
        σ[2] = U[end-1] - U[end];
        IS[2] = σ[2]^2;
        IS[3] = 1/(c[3]^2)*(4/3.*σ[1]^2 + (c[2] - 11/3)*σ[2]*σ[1] +
            (10/3 - 3*c[2] + c[2]^2)*σ[2]^2);
    else # interior
        σ[1] = U[j] - U[j-1];
        σ[2] = U[j+1] - U[j];
        IS[1] = σ[1]^2;
        IS[2] = 13/(12*c[2]^2)*(U[j+1] - 2*U[j] + U[j-1])^2 + 0.25*(U[j+1] - U[j-1])^2;
        IS[3] = σ[2]^2;
    end
end
