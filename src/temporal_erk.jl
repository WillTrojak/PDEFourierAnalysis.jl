abstract type AbstractERKTemporal <: AbstractRKTemporal end

struct ERKTemporal <: AbstractERKTemporal
    A::Array{AbstractFloat, 2}
    b::Vector{AbstractFloat}
    c::Vector{AbstractFloat}

    name::String
    explicit::Bool

    function ERKTemporal(scheme::String)
        A, b, c = get_RK_scheme(scheme)
        explicit = true
        new(A, b, c, scheme, explicit)
    end
end
