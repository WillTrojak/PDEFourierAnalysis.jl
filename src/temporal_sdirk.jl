abstract type AbstractSDIRKTemporal <: AbstractRKTemporal end

mutable struct SDIRKTemporal <: AbstractSDIRKTemporal
    A::Array{AbstractFloat, 2}
    b::Vector{AbstractFloat}
    c::Vector{AbstractFloat}

    i::Integer

    name::String
    explicit::Bool

    function SDIRKTemporal(scheme::String, i::Integer)
        A, b, c = get_RK_scheme(scheme)
        explicit = false
        new(A, b, c, i, scheme, explicit)
    end
end

function pseudo_source(λ::N, dt, X::T) where {T<:AbstractSDIRKTemporal,N<:Number}
    s = 1
    for j=1:X.i-1
        s += X.dt*X.A[X.i,j]*λ*exp(λ*X.c[j]*dt)
    end
    s/(dt*X.A[X.i,X.i])
end

function pseudo_source(Λ::Matrix{N}, dt, X::T) where {T<:AbstractSDIRKTemporal,N<:Number}
    s = zeros(size(Λ)) + UniformScaling(1)
    for j=1:X.i-1
        s += X.dt*X.A[X.i,j]*Λ*exp(Λ*X.c[j]*dt)
    end
    s/(dt*X.A[X.i,X.i])
end

ω(X::T) where {T<:AbstractSDIRKTemporal} = X.A[X.i,X.i]
