abstract type AbstractDualTemporal <: AbstractTemporal end

mutable struct DualTemporal <: AbstractDualTemporal
    Phys::AbstractTemporal
    Pseudo::AbstractTemporal

    m::Integer
    dt::Real
    dτ::Real
    switch
    function DualTemporal(Phys::T, Pseudo::Tau, dt::Real, dτ::Real, pseudo_stage::Integer,
        switch=UniformScaling(1)) where {T<:AbstractTemporal,Tau<:AbstractTemporal}

        if Phys.explicit
            error("Dual time analysis: Physical time scheme not implicit")
        end
        if ~Pseudo.explicit
            error("Dual time analysis: Pseudo time scheme not explicit")
        end

        new(Phys, Pseudo, pseudo_stage, dt, dτ, switch)
    end
end

function amplification_factor(λ, X::T) where {T<:AbstractDualTemporal}
    P = pseudo_amplification_factor(λ, ω(X.Phys), X.dt, X.dτ, X.Pseudo)
    C = pseudo_source_factor(λ, ω(X.Phys), X.dτ, X.Pseudo)
    S = X.switch*pseudo_source(λ, X.dt, X.Phys)
    P^X.m + inv(UniformScaling(1) - P)*(UniformScaling(1) - P^X.m)*C*S
end
