"""
Structure contenant les données de la barre en traction-compression considérée homogène et isotrope

# Paramètres du problème
* L: Longueur [m]
* S: Aire de la section [m²]
* E: Module d'Young [Pa]
* ρ: Masse volumique [kg/m³]

# Paramètres du modèle
* L : Longueur [m]
* m : Masse linéique [kg/m]
* D : Coefficient de raideur [Pa]
"""
@with_kw struct Bar
    L::Float64
    m::Float64
    D::Float64

    function Bar(L::Float64, S::Float64, E::Float64, ρ::Float64)
        m = ρ*S
        D = E*S

        return new(L, m, D)
    end
end

"""
Structure contenant les données de la barre en torsion considérée homogène et isotrope

# Paramètres du constructeur
* L: Longueur [m]
* I: Moment d'inertie [m⁴]
* J: Moment de torsion [m⁴]
* G: Module de Coulomb [Pa]
* ρ: Masse volumique [kg/m³]

# Paramètre de la structure
* L : Longueur [m]
* m : Coefficient d'inertie [kg.m]
* D : Coefficient de raideur [Pa.m²]
"""
@with_kw struct Rod
    L :: Float64
    m::Float64
    D::Float64

    function Rod(L::Float64, I::Float64, J::Float64, G::Float64, ρ::Float64)
        m = ρ*I
        D = G*J

        return new(L, m, D)
    end
end

"""
    eigval(b::Bar, fmax, bc)
    eigval(b::Rod, fmax, bc)

Calcul les fréquences propres d'une barre en traction-compression ou en torsion jusqu'à fmax

# Paramètres
* p: Structure contenant les données relative à la barre
* fₘₐₓ: Fréquence maximale de calcul des déformées modales [Hz]
* bc: Conditions aux limites
    * :CC : Encastrée - Encastrée
    * :CF : Encastrée - Libre
    * :FF : Libre - Libre

# Sorties
* ωₙ: Pulsations propres calculées jusqu'à ωmax = 2π*fmax [Hz]
* kₙ: Vecteur des nombres d'onde modaux
"""
function eigval(b, fₘₐₓ, bc = :CC)
    (; L, m, D) = b

    c = sqrt(D/m)
    ωₘₐₓ = 2π*fₘₐₓ

    ωₙ = Float64[]
    kₙ = Float64[]
    if bc == :CC
        n = 1
        kᵢ = n*π/L
        ωᵢ = c*kᵢ
        while ωᵢ ≤ ωₘₐₓ
            push!(ωₙ, ωᵢ)
            push!(kₙ, kᵢ)
            n += 1
            kᵢ = n*π/L
            ωᵢ = c*kᵢ
        end
    elseif bc === :CF
        n = 1
        kᵢ = (2n - 1)π/L
        ωᵢ = c*kᵢ
        while ωᵢ ≤ ωₘₐₓ
            push!(ωₙ, ωᵢ)
            push!(kₙ, kᵢ)
            n += 1
            kᵢ = (2n - 1)π/L
            ωᵢ = c*kᵢ
        end
    elseif bc == :FF
        n = 0
        kᵢ = 0.
        ωᵢ = 0.
        while ωᵢ ≤ ωₘₐₓ
            push!(ωₙ, ωᵢ)
            push!(kₙ, kᵢ)
            n += 1
            kᵢ = n*π/L
            ωᵢ = c*kᵢ
        end
    else
        error("Boundary conditions not implemented")
    end

    return ωₙ, kₙ
end

"""
    eigmode(b::Bar, kₙ, x, bc)
    eigmode(b::Rod, kₙ, x, bc)

Calcul les déformées propres d'une barre en traction-compression ou en torsion

# Paramètres
* b: Structure contenant les données relative à la barre
* kₙ: Vecteur des nombres d'ondes modaux
* x: Coordonnées des points de calcul des déformées
* bc: Conditions aux limites
    * :CC : Encastrée - Encastrée
    * :CF : Encastrée - Libre
    * :FF : Libre - Libre

# Sorties
* ϕ: Déformées modales normalisées à la masse
"""
function eigmode(b, kₙ, x, bc = :CC)
    (; L, m) = b

    if !isa(x, Array)
        x = collect(x);
    end

    if bc == :CC
        Mₙ = m*L/2

        return sin.(x*kₙ')./sqrt(Mₙ)
    elseif bc == :CF
        Mₙ = m*L/2

        return sin.(x*kₙ')./sqrt(Mₙ)
    elseif bc == :FF
        Mₙ = m*L.*ones(length(n))./2
        Mₙ[1] *= 2.

        return cos.(x*kₙ')./sqrt.(Mₙ')
    else
        error("Boundary conditions not implemented")
    end
end