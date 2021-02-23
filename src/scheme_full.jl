abstract type AbstractFullScheme <: AbstractScheme end

mutable struct FullScheme <: AbstractFullScheme
    X::AbstractSpatial
    T::AbstractTemporal
    function FullScheme(X::Tx, T::Tt) where {Tx<:AbstractSpatia,Tt<:AbstreactTemporal}
        new(X, T)
    end
end
