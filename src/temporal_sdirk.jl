mutable struct SDIRKTemporalScheme <: RKTemporalScheme
    A::Array{AbstractFloat, 2}
    b::Vector{AbstractFloat}
    c::Vector{AbstractFloat}

    i::Integer

    name::String
    explicit::Bool

    function SDIRKTemporalScheme(scheme::String, i::Integer)
        A, b, c = get_RK_scheme(scheme)
        explicit = false
        new(A, b, c, i, scheme, explicit)
    end
end

function pseudo_source(lambda, dt, X::T) where {T<:SDIRKTemporalScheme}
    s = 1
    for j=1:X.i-1
        s += X.dt*X.A[X.i,j]*lambda*exp(lambda*X.c[j]*dt)
    end
    s/(dt*X.A[X.i,X.i])
end

omega(X::T) where {T<:SDIRKTemporalScheme} = X.A[X.i,X.i]
