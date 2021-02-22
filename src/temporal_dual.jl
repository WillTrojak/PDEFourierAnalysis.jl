mutable struct DualTemporalScheme <: TemporalScheme
    Phys::TemporalScheme
    Pseudo::TemporalScheme

    m::Integer
    dt::Real
    dtau::Real
    function DualTemporalScheme(Phys::T, Pseudo::Tau, dt::Real, dtau::Real, pseudo_stage::Integer) where {T<:TemporalScheme,Tau<:TemporalScheme}

        if Phys.explicit
            error("Dual time analysis: Physical time scheme not implicit")
        end
        if ~Pseudo.explicit
            error("Dual time analysis: Pseudo time scheme not explicit")
        end

        new(Phys, Pseudo, dt, dtau, pseudo_stage)
    end
end

function amplification_factor(lambda, X::T) where {T<:DualTemporalScheme}
    P = pseudo_amplification_factor(lambda, omega(X.Phys), X.dt, X.dtau, X.Pseudo)
    C = pseudo_source_factor(lambda, omega(X.Phys), X.dtau, X.Pseudo)
    S = pseudo_source(lambda, X.dt, X.Phys)
    P^X.m + C*S*(1 - P^X.m)/(1 - P)
end
