"""
Structure contenant les données de la poutre en flexion considérée homogène et isotrope

# Paramètres du problème
* L : Longueur [m]
* S : Aire de la section [m²]
* I : Moment quadratique [m⁴]
* E : Module d'Young [Pa]
* ρ : Masse volumique [kg/m³]

# Paramètres du modèle
* L : Longueur [m]
* M : Masse linéique [kg/m]
* D : Raideur de flexion [N.m²]
"""
@with_kw struct Beam
    L :: Float64
    m::Float64
    D::Float64

    function Beam(L::Float64, S::Float64, I::Float64, E::Float64, ρ::Float64)
        m = ρ*S
        D = E*I

        return new(L, m, D)
    end
end

"""
    eigval(b::Beam, fₘₐₓ, bc)

Calcul les fréquences propres d'une poutre en flexion jusqu'à fmax

# Paramètres
* p: Structure contenant les données relative à la barre
* fₘₐₓ: Fréquence maximale de calcul des déformées modales [Hz]
* bc: Conditions aux limites
    * :SS : Appuyée - Appuyée
    * :CC : Encastrée - Encastrée
    * :SC : Appuyée - Encastrée
    * :CF : Encastrée - Libre
    * :SF : Appuyée - Libre
    * :FF : Libre - Libre

# Sorties
* ωₙ: Pulsations propres calculées jusqu'à ωmax = 2π*fmax [Hz]
* kₙ: Vecteur des nombres d'onde modaux
"""
function eigval(b::Beam, fₘₐₓ, bc = :SS)
    (; L, m, D) = b

    c = sqrt(D/m)
    ωₘₐₓ = 2π*fₘₐₓ

    ωₙ = Float64[]
    kₙ = Float64[]
    if bc == :SS
        n = 1
        kᵢ = n*π/L
        ωᵢ = c*kᵢ^2
        while ωᵢ ≤ ωₘₐₓ
            push!(kₙ, kᵢ)
            push!(ωₙ, ωᵢ)
            n += 1
            kᵢ = n*π/L
            ωᵢ = c*kᵢ^2
        end
    elseif bc == :CC
        append!(kₙ, [4.73, 7.85, 11]./L)
        append!(ωₙ, c.*kₙ.^2)

        n = 4
        kᵢ = (2n + 1)π/2L
        ωᵢ = c*kᵢ^2
        while ωᵢ ≤ ωₘₐₓ
            push!(kₙ, kᵢ)
            push!(ωₙ, ωᵢ)
            n += 1
            kᵢ = (2n + 1)π/2L
            ωᵢ = c*kᵢ.^2
        end
    elseif bc == :CS
        append!(kₙ, [3.92, 7.07, 10.2]./L)
        append!(ωₙ, c.*kₙ.^2)

        n = 4
        kᵢ = (4n + 1)π/4L
        ωᵢ = c*kᵢ^2
        while ωᵢ ≤ ωₘₐₓ
            push!(kₙ, kᵢ)
            push!(ωₙ, ωᵢ)
            n += 1
            kᵢ = (4n + 1)π/4L
            ωᵢ = c*kᵢ^2
        end
    elseif bc == :CF
        append!(kₙ, [1.87,  4.73, 7.85, 11]./L)
        append!(ωₙ, c.*kₙ.^2)

        n = 5
        kᵢ = (2n + 1)π/2L
        ωᵢ = c*kᵢ^2
        while ωᵢ ≤ ωₘₐₓ
            push!(kₙ, kᵢ)
            push!(ωₙ, ωᵢ)
            n += 1
            kᵢ = (2n + 1)π/2L
            ωᵢ = c*kᵢ^2
        end
    elseif bc == :SF
        append!(kₙ, [3.92, 7.07, 10.2]./L)
        append!(ωₙ, c.*kₙ.^2)

        n = 4
        kᵢ = (4n + 1)π/4L
        ωᵢ = c*kᵢ^2
        while ωᵢ ≤ ωₘₐₓ
            push!(kₙ,kᵢ)
            push!(ωₙ, ωᵢ)
            n += 1
            kᵢ = (4n + 1)π/4L
            ωᵢ = c*kᵢ^2
        end
    elseif bc == :FF
        append!(kₙ, [0., 0., 4.73, 7.85, 11]./L)
        append!(ωₙ, c.*kₙ.^2)

        n = 4
        kᵢ = (2n + 1)π/2L
        ωᵢ = c*kᵢ^2
        while ωᵢ ≤ ωₘₐₓ
            push!(kₙ, kᵢ)
            push!(ωₙ, ωᵢ)
            n += 1
            kᵢ = (2n + 1)π/2L
            ωᵢ = c*kᵢ^2
        end
    else
        error("Boundary conditions not implemented")
    end

    return ωₙ, kₙ
end

"""
    eigmode(b::Beam, kₙ, x, bc)

Calcul les déformées propres d'une poutre en flexion

# Paramètres
* b: Structure contenant les données relative à la poutre
* kₙ: Vecteur des nombres d'ondes modaux
* x: Coordonnées des points de calcul des déformées
* bc: Conditions aux limites
    * :SS : Appuyée - Appuyée
    * :CC : Encastrée - Encastrée
    * :CS : Encastrée - Appuyée
    * :CF : Encastrée - Libre
    * :SF : Appuyée - Libre
    * :FF : Libre - Libre

# Sorties
* ϕ: Déformées modales normalisées à la masse
"""
function eigmode(b::Beam, kₙ, x, bc = :SS)
    (; L, m) = b

    if !isa(x, Array)
        x = collect(x);
    end

    if bc == :SS
        Mₙ = m*L/2.
        ϕₙ = sin.(x*kₙ')
    elseif bc == :CC
        Mₙ = @. m*((-kₙ*L*cos(2kₙ*L) + kₙ*L*cosh(2kₙ*L) + cosh(kₙ*L)^2*sin(2kₙ*L) + 2cos(kₙ*L)*sinh(kₙ*L) - 2sin(kₙ*L)*(cosh(kₙ*L) + 2kₙ*L*sinh(kₙ*L)) - cos(kₙ*L)^2*sinh(2kₙ*L))/(2kₙ*(cos(kₙ*L) - cosh(kₙ*L))^2*(sin(kₙ*L) - sinh(kₙ*L))^2))

        ϕₙ = @. ((cosh(x*kₙ') - cos(x*kₙ'))/(cosh(kₙ'*L) - cos(kₙ'*L))) - ((sinh(x*kₙ') - sin(x*kₙ'))/(sinh(kₙ'*L) - sin(kₙ'*L)))
    elseif bc == :CS
        Mₙ = @. m*((-kₙ*L*cos(2kₙ*L) + kₙ*L*cosh(2kₙ*L) + cosh(kₙ*L)^2*sin(2kₙ*L) + 2cos(kₙ*L)*sinh(kₙ*L) - 2sin(kₙ*L)*(cosh(kₙ*L) + 2kₙ*L*sinh(kₙ*L)) - cos(kₙ*L)^2*sinh(2kₙ*L))/(2kₙ*(cos(kₙ*L) - cosh(kₙ*L))^2*(sin(kₙ*L) - sinh(kₙ*L))^2))

        ϕₙ = @. ((cosh(x*kₙ') - cos(x*kₙ'))/(cosh(kₙ'*L) - cos(kₙ'*L))) - ((sinh(x*kₙ') - sin(x*kₙ'))/(sinh(kₙ'*L) - sin(kₙ'*L)))

    elseif bc == :CF
        Mₙ = @. m*((-2kₙ*L*cos(2kₙ*L) - 7sin(2kₙ*L) + cosh(2kₙ*L)*(2kₙ*L + 3sin(2kₙ*L)) - 2cosh(kₙ*L)*(3sin(kₙL) + sin(3kₙ*L)) - sin(4kₙ*L) + 2(3(cos(kₙ*L) + cos(3kₙ*L)) - 4kₙ*L*sin(kₙ*L))*sinh(kₙ*L) + 6cos(kₙ*L)^2*sinh(2kₙ*L))/(4kₙ*(cos(kₙ*L) + cosh(kₙ*L))^2*(sin(kₙ*L) - sinh(kₙ*L)).^2))

        ϕₙ = @. ((cosh(x*kₙ') - cos(x*kₙ'))/(cosh(kₙ'*L) + cos(kₙ'*L))) - ((sinh(x*kₙ') - sin(x*kₙ'))/(sinh(kₙ'*L) + sin(kₙ'*L)))
    elseif bc == :SF
        Mₙ = @. m*((-3/tan(kₙ*L) + 3/tanh(kₙ*L) + kₙ*L/(sin(kₙ*L)^2) - kₙ*L/(sinh(kₙ*L)^2))/2kₙ);

        ϕₙ = @. (sin(x*kₙ')/sin(kₙ'*L)) + (sinh(x*kₙ')/sinh(kₙ'*L))
    elseif bc == :FF
        Mₙ = @. m*((-kₙ*L*cos(2kₙ*L) + kₙ*L*cosh(2kₙ*L) + 6cosh(kₙ*L)*sin(kₙ*L) - 3cosh(kₙ*L)^2*sin(2kₙ*L) - 2(3cos(kₙ*L) + 2kₙ*L*sin(kₙ*L))*sinh(kₙ*L) + 3cos(kₙ*L)^2*sinh(2kₙ*L))/(2kₙ*(cos(kₙ*L) - cosh(kₙ*L))^2*(sin(kₙ*L) - sinh(kₙ*L))^2))
        Mₙ[1:2] .= m*L

        ϕₙ = @. ((cosh(x*kₙ') + cos(x*kₙ'))/(cosh(kₙ'*L) - cos(kₙ'*L))) - ((sinh(x*kₙ') + sin(x*kₙ'))/(sinh(kₙ'*L) - sin(kₙ'*L)))

        ϕₙ[:, 1] .= 1.
        ϕₙ[:, 2] = x .- L/2
    else
        error("Boundary conditions not implemented")
    end

    return @. ϕₙ/sqrt(Mₙ')
end