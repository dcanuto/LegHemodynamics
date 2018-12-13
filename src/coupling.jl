function coupling!(system::LegSystem,n::Int64,h::Float64,terms::Vector{Int64},uprox::Float64)
    # println("Distal coupling:")
    CVModule.coupledistal!(system,n,h,terms);
    # println("Proximal coupling:")
    CVModule.coupleproximal!(system,n,h,uprox);
end
