"""
Structure containing data for the time solver

# Parameters
* K: Stiffness matrix
* M: Mass matrix
* C: Damping matrix
* t: Collection of calculation time steps
"""
struct LinearTimeProblem
    K
    M
    C
    F
    t
end

"""
Structure containing problem solutions

# Parameters
* D: Displacement matrix
* V: Velocity matrix
* A: Acceleration matrix
"""
struct TimeSolution
    D :: Matrix{Float64}
    V :: Matrix{Float64}
    A :: Matrix{Float64}
end


# Solvers
struct CentralDiff end # Différences finies centrées

struct RK4 end # Runge-Kutta d'ordre 4

@with_kw struct FoxGoodwin
    αf::Float64 = 0.
    αₘ::Float64 = 0.
    γ₀::Float64 = 0.5
    β₀::Float64 = 1/12
end # Fox-Goodwin

@with_kw struct LinearAcceleration
    αf::Float64 = 0.
    αₘ::Float64 = 0.
    γ₀::Float64 = 0.5
    β₀::Float64 = 1/6
end # Linear Acceleration

struct Newmark
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64

    Newmark(; γ₀ = 0.5, β₀ = 0.25) = new(0., 0., γ₀, β₀)
end # Newmark

struct HHT
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64

    function HHT(; γ₀ = 0.5, β₀ = 0.25, ρ = 1., αf = Inf)
        if αf ≠ Inf
            (0. < αf < 1/3) ? error("αf must be in [0, 1/3[") : nothing
            return new(αf, 0., γ₀, β₀)
        else
            ρ < 0.5 ? error("ρ must be in [0.5, 1]") : nothing
            return new((1. - ρ)/(1. + ρ), 0., γ₀, β₀)
        end
    end
end # Hilber-Hughes-Taylor

struct WBZ
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64

    function WBZ(; γ₀ = 0.5, β₀ = 0.25, ρ = 1., αₘ = Inf)
        if αₘ ≠ Inf
            (αₘ ≤ 0.5) ? error("αₘ must be ≤ 0.5") : nothing
            return new(0., αₘ, γ₀, β₀)
        else
            (ρ > 1.) ? error("ρ must be in [0, 1]") : nothing
            return new(0., (ρ - 1.)/(ρ + 1.), γ₀, β₀)
        end
    end
end # Wood-Bossak-Zienkiewicz

struct GeneralizedAlpha
    αf::Float64
    αₘ::Float64
    γ₀::Float64
    β₀::Float64

    function GeneralizedAlpha(; γ₀ = 0.5, β₀ = 0.25, ρ = 1., αf = Inf, αₘ = Inf)
        if αf ≠ Inf && αₘ ≠ Inf
            (αₘ ≤ αf ≤ 0.5) ? error("αₘ ≤ αf ≤ 0.5") : nothing
        else
            (ρ > 1.) ? error("ρ must be in [0, 1]") : nothing

            αf = ρ/(ρ + 1.)
            αₘ = 3αf - 1.
        end

        return new(αf, αₘ, γ₀, β₀)
    end
end # Generalized-α

@with_kw struct MidPoint
    αf::Float64 = 0.5
    αₘ::Float64 = 0.5
    γ₀::Float64 = 0.5
    β₀::Float64 = 0.25
end # Mid-Point rule

function solve(prob::LinearTimeProblem, u0, alg = GeneralizedAlpha())
    if alg isa CentralDiff # Central-difference
        return _CentralDiff(prob, u0)
    elseif alg isa RK4 # Fourth-order Runge-Kutta
        return _RK4(prob, u0)
    else
        return _NK(prob, u0, alg)
    end
end

# Central-difference algorithm
function _CentralDiff(prob::LinearTimeProblem, u0)
    (; K, M, C, F, t) = prob

    nt = length(t)
    h = (maximum(t) - minimum(t))/(nt - 1) # Pas de temps
    Nddl = size(K)[1]

    # Initialisation des matrices de résultat
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Calcul de l'accélération initiale
    D[:, 1] = u0.D₀
    V[:, 1] = u0.V₀

    rhs0 = F[:, 1] - C*V[:, 1] - K*D[:, 1]
    LU = lu(M)
    A[:, 1] = LU\rhs0

    D_1 = D[:, 1] - h.*V[:, 1] + (h^2 .*A[:, 1]./4)
    p = Progress(nt - 1; desc = "Central difference...", color = :black, barlen = 75, showspeed = true)
    for n in 1:nt-1
        next!(p)
        if n == 1
            D[:, n+1] = 2D[:, n] - D_1 + (h^2 .*A[:, n])
        else
            D[:, n+1] = 2D[:, n] - D[:, n-1] + (h^2 .*A[:, n])
        end

        V[:, n+1] = (D[:, n+1] - D[:, n])./h

        rhs = F[:, n+1] - C*V[:, n+1] - K*D[:, n+1]
        A[:, n+1] = LU\rhs
    end

    return TimeSolution(D, V, A)
end

@views function _RK4(prob::LinearTimeProblem, u0)
    (; K, M, C, F, t) = prob

    nt = length(t)
    h = (maximum(t) - minimum(t))/(nt - 1) # Pas de temps
    Nddl = size(K)[1]

    # Initialisation des matrices de résultat
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Calcul de l'accélération initiale
    D[:, 1] = u0.D₀
    V[:, 1] = u0.V₀

    rhs0 = F[:, 1] - C*V[:, 1] - K*D[:, 1]
    LU = lu(M)
    A[:, 1] = LU\rhs0

    p = Progress(nt - 1; desc = "RK4...", color=:black, barlen=75, showspeed=true)
    for n = 1:nt-1
        next!(p)
        Fn_2 = (F[:, n+1] + F[:, n])./2

        k₁ = A[:, n]

        CK = (C + h.*K./2)*V[:, n]
        KD = K*D[:, n]
        k₂ = Fn_2 - CK - h.*(C*k₁)./2 - KD
        k₃ = Fn_2 - CK - h.*(C*k₂)./2 - (h^2 .*(K*k₁)./4) - KD
        k₄ = F[:, n+1] - (C + h.*K)*V[:, n] - h.*(C*k₃) - (h^2 .*(K*k₂)./2) - KD

        D[:, n+1] = D[:, n] + h.*V[:, n] + (h^2 .*(k₁ + (LU\(k₂ + k₃)))./6)
        V[:, n+1] = V[:, n] + h.*(k₁ + (LU\(2k₂ + 2k₃ + k₄)))./6
        A[:, n+1] = LU\(F[:, n+1] - C*V[:, n+1] - K*D[:, n+1])
    end

    return TimeSolution(D, V, A)
end

function _NK(prob::LinearTimeProblem, u0, alg) # Newmark-Family
    (; K, M, C, F, t) = prob
    (; αf, αₘ, γ₀, β₀) = alg

    if alg isa FoxGoodwin
        desc = "Fox-Goodwin..."
    elseif alg isa LinearAcceleration
        desc = "Linear acceleration..."
    elseif alg isa Newmark
        desc = "Newmark..."
    elseif alg isa HHT
        desc = "HHT..."
    elseif alg isa WBZ
        desc = "WBZ..."
    elseif alg isa GeneralizedAlpha
        desc = "Generalized-α..."
    elseif alg isa MidPoint
        desc = "Mid-point rule..."
    end

    nt = length(t)
    h = (maximum(t) - minimum(t))/(nt - 1) # Pas de temps
    Nddl = size(K)[1]

    # Initialisation des matrices de résultat
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Calcul de l'accélération initiale
    D[:, 1] = u0.D₀
    V[:, 1] = u0.V₀

    rhs0 = F[:, 1] - C*V[:, 1] - K*D[:, 1]
    A[:, 1] = M\rhs0

    # Paramètre du schéma de Newmark
    γ = γ₀*(1. + 2αf - 2αₘ)
    β = β₀*(1. + αf - αₘ)^2

    a₁ = (1. - γ)*h
    a₂ = (0.5 - β)*h^2
    a₃ = γ*h
    a₄ = β*h^2

    b₈ = 1. - αf
    b₁ = b₈*a₁
    b₂ = b₈*a₂
    b₃ = b₈*a₃
    b₄ = b₈*a₄
    b₅ = b₈*h
    b₆ = 1. - αₘ
    b₇ = αₘ
    b₉ = αf

    # Calcul de la raideur effective
    S = @. b₆*M + b₃*C + b₄*K
    LU = lu(S)

    p = Progress(nt - 1; desc = desc, color = :black, barlen = 75, showspeed = true)
    for n = 1:nt-1
        next!(p)
        rhs = b₈.*F[:, n+1] + b₉.*F[:, n] - C*(b₁.*A[:, n] + V[:, n]) - K*(b₂.*A[:, n] + b₅.*V[:, n] + D[:, n]) - b₇.*(M*A[:, n])

        A[:, n+1] = LU\rhs
        V[:, n+1] = @. V[:, n] + a₁*A[:, n] + a₃*A[:, n+1]
        D[:, n+1] = @. D[:, n] + h*V[:, n] + a₂*A[:, n] + a₄*A[:, n+1]
    end

    return TimeSolution(D, V, A)
end