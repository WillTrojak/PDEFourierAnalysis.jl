mutable struct FullScheme <: AbstractScheme
    X::AbstractSpatial
    T::abstractTemporal
    function FullScheme(X::Tx, T::Tt) where {Tx<:AbstractSpatia,Tt<:AbstreactTemporal}
        new(X, T)
    end
end
