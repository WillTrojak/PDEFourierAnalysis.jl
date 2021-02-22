struct ERKTemporalScheme <: RKTemporal
    A::Array{AbstractFloat, 2}
    b::Vector{AbstractFloat}
    c::Vector{AbstractFloat}

    name::String
    explicit::Bool

    function ERKTemporalScheme(scheme::String)
        A, b, c = get_RK_scheme(scheme)
        explicit = true
        new(A, b, c, scheme, explicit)
    end
end
