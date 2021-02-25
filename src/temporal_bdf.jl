abstract type AbstractBDFTemporal <: AbstractTemporal end

mutable struct BDFTemporal <: AbstractBDFTemporal
    B::Vector{Number}

    name::String
    explicit::Bool
    function BDFTemporal(scheme::String)
        B = @eval $(Meta.parse(string("B_", scheme, "()")))
        explicit = false
        new(vec(B), scheme, explicit)
    end
end

function pseudo_source(λ::N, dt, X::T) where {T<:AbstractSDIRKTemporal,N<:Number}
    s = 0
    for j=0:length(X.B)-1
        s -= X.B[j+2]*exp(-j*λ*dt)
    end
    s/(dt*X.B[1])
end

function pseudo_source(Λ::Matrix{N}, dt, X::T) where {T<:AbstractSDIRKTemporal,N<:Number}
    s = zeros(size(Λ))
    for j=0:length(X.B)-1
        s -= X.B[j+2]*exp(-j*Λ*dt)
    end
    s/(dt*X.B[1])
end

ω(X::T) where {T<:AbstractBDFTemporal} = X.B[1]
