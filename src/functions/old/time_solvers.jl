"""
Structure contenant les données alimentant le solver temporel

# Paramètres
* K: Matrice de raideur
* M: Matrice de masse
* C: Matrice d'amortissement
* t: Ensemble des pas de temps de calcul
"""
struct TimeProblem
    K
    M
    C
    F
    t
end

"""
Structure contenant les solutions du problème

# Paramètres
* D: Matrice d'évolution des déplacements
* V: Matrice d'évolution des vitesses
* A: Matrice d'évolution des accélérations
"""
struct TimeSolution
    D :: Matrix{Float64}
    V :: Matrix{Float64}
    A :: Matrix{Float64}
end

"""
    solve(prob::TimeProblem, CI, type = :CD, ρ = 1.)

Résolution du problème par schéma d'intégration temporelle

# Paramètres
* prob : Structure contenant les données du problème
* CI : NamedTuples des conditions initiales
    - D₀ : Déplacement initial - Dimension Nddl
    - V₀ : Vitesse initiale - Dimension Nddl
* type : Méthode d'intégration
* ρ : Rayon spectral (nécessaire pour les méthodes HHT, WBZ et Gα)
* form : Méthode de calcul pour les schémas de la famille Newmark (:ST - standard, :AB - méthode d'Arnorl and Brüels)

# Méthodes implémentées
* :CD : Différences finies centrées
* :RK4 : Runge-Kutta d'ordre 4
* :FG : Fox-Goodwin
* :LA : Linear Acceleration
* :NK : Newmark
* :MP : Mid-Point rule
* :HHT : Hilber-Hughes-Taylor
* :WBZ : Wood-Bossak-Zienkiewicz
* :Gα : Generalized-alpha

# Sortie
* sol : Structure de données contenant les matrices de déplacement, vitesse et accélération - Dimension Nddl x nt
"""
function solve(prob::TimeProblem, CI, type = :CD; ρ = 1., form = :ST)
    keys = (:CD, :RK4, :FG, :LA, :NK, :MP, :HHT, :WBZ, :Gα)

    if !(type in keys)
        error("Method not implemented")
    end

    if type == :CD
        sol = central_diff(prob, CI)
    elseif type == :RK4
        sol = RK4(prob, CI)
    elseif type in keys[3:end]
        if form == :ST
            sol = newmark_family(prob, CI, type, ρ)
        else
            sol = newmark_family_AB(prob, CI, type, ρ)
        end
    end

    return sol
end

function central_diff(prob::TimeProblem, CI)
    (; K, M, C, F, t) = prob

    nt = length(t)
    h = (maximum(t) - minimum(t))/(nt - 1) # Pas de temps
    Nddl = size(K)[1]

    # Initialisation des matrices de résultat
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Calcul de l'accélération initiale
    D[:, 1] = CI.D₀
    V[:, 1] = CI.V₀

    rhs0 = F[:, 1] - C*V[:, 1] - K*D[:, 1]
    LU = lu(M)
    A[:, 1] = LU\rhs0

    D_1 = D[:, 1] - h.*V[:, 1] + (h^2 .*A[:, 1]./4)
    p = Progress(nt - 1; color = :black, barlen = 75, showspeed = true)
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

@views function RK4(prob::TimeProblem, CI)
    (; K, M, C, F, t) = prob

    nt = length(t)
    h = (maximum(t) - minimum(t))/(nt - 1) # Pas de temps
    Nddl = size(K)[1]

    # Initialisation des matrices de résultat
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Calcul de l'accélération initiale
    D[:, 1] = CI.D₀
    V[:, 1] = CI.V₀

    rhs0 = F[:, 1] - C*V[:, 1] - K*D[:, 1]
    LU = lu(M)
    A[:, 1] = LU\rhs0

    p = Progress(nt - 1; color=:black, barlen=75, showspeed=true)
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

function newmark_family(prob::TimeProblem, CI, type, ρ)
    (; K, M, C, F, t) = prob

    nt = length(t)
    h = (maximum(t) - minimum(t))/(nt - 1) # Pas de temps
    Nddl = size(K)[1]

    # Initialisation des matrices de résultat
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Calcul de l'accélération initiale
    D[:, 1] = CI.D₀
    V[:, 1] = CI.V₀

    rhs0 = F[:, 1] - C*V[:, 1] - K*D[:, 1]
    A[:, 1] = M\rhs0

    # Choix de la méthode
    if type == :FG
        αf = 0.
        αₘ = 0.
        γ₀ = 0.5
        β₀ = 1/12
    elseif type == :LA
        αf = 0.
        αₘ = 0.
        γ₀ = 0.5
        β₀ = 1/6
    else
        γ₀ = 0.5
        β₀ = 0.25
        if type == :NK
            αf = 0.
            αₘ = 0.
        elseif type == :HHT
            if ρ < 0.5
                error("ρ must be in [0.5, 1]")
            end
            αf = (1. - ρ)/(1. + ρ)
            αₘ = 0.

            if (0. < αf < 1/3)
                error("αf must be in [, 1/3[")
            end
        elseif type == :WBZ
            αf = 0.
            if ρ > 1.
                error("ρ ∈ [0, 1]")
            end
            αₘ = (ρ - 1.)/(ρ + 1.)

            if αₘ > 0.5
                error("αₘ must be ≤ 0.5")
            end
        elseif type == :Gα
            αf = ρ/(ρ + 1.)
            αₘ = 3αf - 1.

            if (αₘ > 0.5) || (αf > 0.5)
                error("αₘ ≤ αf ≤ 0.5")
            end
        elseif type == :MP
            αf = 0.5
            αₘ = 0.5
        end
    end

    # Paramètre du schéma de Newmark
    γ = γ₀*(1. +2αf - 2αₘ)
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

    p = Progress(nt - 1; color = :black, barlen = 75, showspeed = true)
    for n = 1:nt-1
        next!(p)
        rhs = b₈.*F[:, n+1] + b₉.*F[:, n] - C*(b₁.*A[:, n] + V[:, n]) - K*(b₂.*A[:, n] + b₅.*V[:, n] + D[:, n]) - b₇.*(M*A[:, n])

        A[:, n+1] = LU\rhs
        V[:, n+1] = @. V[:, n] + a₁*A[:, n] + a₃*A[:, n+1]
        D[:, n+1] = @. D[:, n] + h*V[:, n] + a₂*A[:, n] + a₄*A[:, n+1]
    end

    return TimeSolution(D, V, A)
end

@views function newmark_family_AB(prob::TimeProblem, CI, type, ρ)
    # Implémentation alternative due à Arnold and Brüls (2007)
    (; K, M, C, F, t) = prob

    nt = length(t)
    h = (maximum(t) - minimum(t))/(nt - 1) # Pas de temps
    Nddl = size(K)[1]

    # Initialisation des matrices de résultat
    D = zeros(Nddl, nt)
    V = zeros(Nddl, nt)
    A = zeros(Nddl, nt)

    # Calcul de l'accélération initiale
    D[:, 1] = CI.D₀
    V[:, 1] = CI.V₀

    rhs0 = F[:, 1] - C*V[:, 1] - K*D[:, 1]
    A[:, 1] = M\rhs0

    # Choix de la méthode
    if type == :FG
        αf = 0.
        αₘ = 0.
        γ₀ = 0.5
        β₀ = 1/12
    elseif type == :LA
        αf = 0.
        αₘ = 0.
        γ₀ = 0.5
        β₀ = 1/6
    else
        γ₀ = 0.5
        β₀ = 0.25
        if type == :NK
            αf = 0.
            αₘ = 0.
        elseif type == :HHT
            if ρ < 0.5
                error("ρ must be in [0.5, 1]")
            end
            αf = (1. - ρ)/(1. + ρ)
            αₘ = 0.

            if (0. < αf < 1/3)
                error("αf must be in [0, 1/3[")
            end
        elseif type == :WBZ
            αf = 0.
            if ρ > 1.
                error("ρ ∈ [0, 1]")
            end
            αₘ = (ρ - 1.)/(ρ + 1.)

            if αₘ > 0.5
                error("αₘ must be ≤ 0.5")
            end
        elseif type == :Gα
            αf = ρ/(ρ + 1.)
            αₘ = 3αf - 1.

            if (αₘ > 0.5) || (αf > 0.5)
                error("αₘ ≤ αf ≤ 0.5")
            end
        elseif type == :MP
            αf = 0.5
            αₘ = 0.5
        end
    end

    # Paramètre du schéma de Newmark
    γ = γ₀*(1. +2αf - 2αₘ)
    β = β₀*(1. + αf - αₘ)^2

    a₁ = 1/β/h^2
    a₂ = 1/β/h
    a₃ = 1. - 1/2β

    b₀ = γ/β
    b₁ = γ*a₂
    b₂ = 1. - b₀
    b₃ = h*(1. - b₀/2)

    c₀ = (1. - αₘ)/(1. - αf)
    c₁ = a₁*c₀
    c₂ = a₂*c₀
    c₃ = αₘ/(1. - αₘ) + a₃*c₀
    c₄ = αf/(αf - 1.)

    # Calcul de la raideur effective
    S = c₁.*M + b₁.*C + K
    LU = lu(S)

    qₙ = A[:, 1]
    p = Progress(nt - 1; color=:black, barlen=75, showspeed=true)
    for n = 1:nt-1
        next!(p)
        rhs = F[:, n+1] - M*(c₄.*A[:, n] + c₃.*qₙ - c₂.*V[:, n] - c₁.*D[:,n]) - C*(b₃.*qₙ + b₂.*V[:, n] - b₁.*D[:, n])

        D[:, n+1] = LU\rhs

        Dtemp = D[:, n+1] - D[:, n]
        V[:, n+1] = @. b₁*Dtemp + b₂*V[:, n] + b₃*qₙ
        A[:, n+1] = @. c₁*Dtemp - c₂*V[:, n] + c₃*qₙ + c₄*A[:, n]

        qₙ = @. a₁*Dtemp - a₂*V[:, n] + a₃*qₙ
    end

    return TimeSolution(D, V, A)
end