"""
Structure contenant les données de la plaque en flexion considérée homogène et isotrope

# Paramètres du problème
* L : Longueur [m]
* b : Largeur [m]
* E : Module d'Young [Pa]
* ρ : Masse volumique [kg/m³]
* ν : Coefficient de Poisson

# Paramètres du modèle
* L : Longueur [m]
* b : Largeur [m]
* m : Masse surfacique [kg/m²]
* D : Raideur de flexion [N.m]
"""
@with_kw struct Plate
    L::Float64
    b::Float64
    m::Float64
    D::Float64

    function Plate(L::Float64, b::Float64, h::Float64, E::Float64, ρ::Float64, ν::Float64)

        m = ρ*h
        D = E*h^3/12/(1. - ν^2.)

        return new(L, b, m, D)
    end
end

"""
    eigval(p::Plate, fₘₐₓ)

Calcul les fréquences propres d'une plaque simplement appuyée jusqu'à fmax

# Paramètres
    * p : Structure contenant les données relative à la plaque
    * fₘₐₓ : Fréquence maximale de calcul des déformées modales [Hz]

# Sorties
    * ωₘₙ : Pulsations propres calculées jusqu'à ωmax = 2π*fmax [Hz]
    * kₘₙ : Matrice des nombres d'ondes modaux
"""
function eigval(p::Plate, fₘₐₓ)
   (; L, b, m, D) = p

    c = sqrt(D/m)
    ωₘₐₓ = 2π*fₘₐₓ

    m = 1
    n = 1
    kₘ = m*π/L
    kₙ = n*π/b
    ωᵢ = c*(kₘ^2 + kₙ^2)

    ωₘₙ = Float64[]
    kₘₙ = Float64[]
    ind = Int64[]
    # Boucle sur m
    while ωᵢ ≤ ωₘₐₓ
        # Boucle sur n
        while ωᵢ ≤ ωₘₐₓ
            push!(ωₘₙ,  ωᵢ)
            append!(kₘₙ, [kₘ, kₙ])
            append!(ind, [m, n])
            n += 1
            kₘ = m*π/L
            kₙ = n*π/b
            ωᵢ = c*(kₘ^2 + kₙ^2)
        end

        m += 1
        n = 1
        kₘ = m*π/L
        kₙ = n*π/b
        ωᵢ = c*(kₘ^2 + kₙ^2)
    end

    kₘₙ = reshape(kₘₙ, (2, Int(length(kₘₙ)/2)))
    ind = reshape(ind, (2, Int(length(ind)/2)))
    pos = sortperm(ωₘₙ)

    return ωₘₙ[pos], kₘₙ[:, pos], ind[:, pos]
end

"""
    eigmode(p::Plate, kₘₙ, X, Y)

Calcul les déformées propres d'une plaque simplement appuyée

# Paramètres
    * p : Structure contenant les données relative à la plaque
    * kₘₙ : Matrice des nombres d'onde modaux
    * (X, Y): Coordonnées des points de calcul des déformées

# Sorties
    * ϕ: Déformées modales normalisées à la masse
"""
@views function eigmode(p::Plate, kₘₙ, X, Y)
    (; L, b, m) = p
    Mₙ = m*L*b/4

    return sin.(X*kₘₙ[1, :]').*sin.(Y*kₘₙ[2, :]')./sqrt(Mₙ)
end