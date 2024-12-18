"""
Structure contenant les données alimentant le solver modal de calcul d'une FRF

# Paramètres
* ωₙ : Pulsations de résonance
* ξₙ : Amortissements modaux
* ϕₑ : Déformée modale aux points d'excitation
* ϕₒ : Déformée modale aux points d'observation
* freq : Fréquences d'intérêt

# Note
Les déformées modales doivent être normées à la masse
"""
struct ModalFRF
    ωₙ :: Vector{Float64}
    ξₙ
    ϕₑ
    ϕₒ
    freq
end

"""
Structure contenant les données alimentant le solver direct de calcul d'une FRF

# Paramètres
* K : Matrice de raideur
* M : Matrice de masse
* C : Matrice d'amortissement
* exc_dofs : Degrés de liberté d'excitation
* obs_dofs : Degrés de liberté d'observation
* freq : Fréquences d'intérêt
"""
struct DirectFRF
    K
    M
    C
    exc_dofs
    obs_dofs
    freq
end

"""
    frf(m::ModalFRF)
    frf(m::DirectFRF)

Calcul de la matrice des FRF par approche modale ou directe

# Paramètre
* m : Structure contenant les données du problème

# Sortie
* FRF : Matrice des FRF
"""
function frf(m::ModalFRF, type = :dis)
    # Initialisation
    (; ωₙ, ξₙ, ϕₑ, ϕₒ, freq) = m
    Nₑ = size(ϕₑ, 1)
    Nₒ = size(ϕₒ, 1)
    Nf = length(freq)

    FRF = zeros(Complex{Float64}, Nₒ, Nₑ, Nf)

    ωf = 2π*freq
    p = Progress(Nf, color=:black, barlen=75, showspeed=true)
    @inbounds for (f, ω) in enumerate(ωf)
        next!(p)
        M = spdiagm(@. 1/(ωₙ^2 - ω^2 + 2im*ξₙ*ωₙ*ω))
        FRF[:, :, f] = ϕₑ*M*ϕₒ'

        if type == :vel
            FRF[:, :, f] *= 1im*ω
        elseif type == :acc
            FRF[:, :, f] *= -ω^2
        end
    end

    return FRF
end