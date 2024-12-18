"""
    excitation(param, t)

Calcul de différents type de signaux d'excitation

# Paramètres
* param : NamedTuple ou struct contenant les champs suivants :
    * type : type de l'excitation
        1. :triangle
        2. :rectangle
        3. :marteau
        4. :srect (smooth rectangle)
    * F₀ : Amplitude de l'effort [N]
    * td : Instant de déclenchement de l'effort [s]
    * duree : Durée de l'excitation [s]
        *note : Obligatoire pour les type 1, 2 et 4, inutile pour le type 3*
    * k et θ : Paramètre de forme [-] et paramètre d'intensité [s]
        *note : Obligatoire pour le type 3, inutile pour les autres types*
    * tm : Temps de montée entre 0 et F₀ [s]
        *note : Obligatoire pour le type 4, inutile pour les autres*

# Sortie
* F : Vecteur d'évolution de l'excitation en fonction du temps [N]
"""
function excitation(param, t)
    (; type, F₀, td) = param

    if type == :triangle
        (; duree) = param
        F = triangle(F₀, td, duree, t)
    elseif type == :rectangle
        (; duree) = param
        F = rectangle(F₀, td, duree, t)
    elseif type == :marteau
        (; k, θ) = param
        F = marteau(F₀, td, k, θ, t)
    elseif type == :srect
        (; tm, duree) = param
        F = smooth_rect(F₀, td, tm, duree, t)
    else
        error("Excitation type not implemented")
    end
end

"""
    triangle(F₀, td, duree, t)

Calcul du vecteur d'excitation pour un signal triangulaire

# Paramètres
* F0 : Amplitude de l'excitation [N]
* td : Instant de déclenchement de l'effort [s]
* duree : Durée de l'excitation [s]
* t : Discrétisation temporelle [s]

# Sortie
* Ft : Vecteur d'évolution de l'excitation en fonction du temps [N]
"""
function triangle(F₀, td, duree, t)
    Ft = zeros(length(t))

    tm = (2td + duree)/2.
    pos_deb = argmin((t .- td).^2.)
    pos_milieu = argmin((t .- tm).^2.)
    pos_fin = argmin((t .- td .- duree).^2.)
    amp = 2F₀/duree

    Ft[pos_deb:pos_milieu] = amp*(t[pos_deb:pos_milieu] .- td)
    Ft[pos_milieu + 1:pos_fin] = F₀ .- amp*(t[pos_milieu + 1:pos_fin] .- tm)

    return Ft
end

"""
    rectangle(F₀, td, duree, t)

Calcul du vecteur d'excitation pour un signal de type créneau

# Paramètres
* F0 : Amplitude de l'excitation [N]
* td : Instant de déclenchement de l'effort [s]
* duree : Durée de l'excitation [s]
* t : Discrétisation temporelle [s]

# Sortie
* Ft : Vecteur d'évolution de l'excitation en fonction du temps [N]
"""
function rectangle(F₀, td, duree, t)
    Ft = zeros(length(t))

    pos_deb = argmin((t .- td).^2.)
    pos_fin = argmin((t .- td .- duree).^2.)

    pos_exc_t = findall(t[pos_deb] .≤ t .≤ t[pos_fin])

    Ft[pos_exc_t] .= F₀

    return Ft
end

"""
    marteau(F₀, td, k, θ, t)

Calcul du vecteur d'excitation pour un signal de type coup de marteau
L'excitation par coup de marteau est supposée avoir la forme d'une loi gamma
de paramètres (k, theta)

# Paramètres
* F0 : Amplitude de l'excitation [N]
* td : Instant de déclenchement de l'effort [s]
* k : Paramètre de forme
* θ : Paramètre d'intensité [s]
* t : Discrétisation temporelle [s]

# Sortie
* Ft : Vecteur d'évolution de l'excitation en fonction du temps [N]
"""
function marteau(F₀, td, k, θ, t)
    Ft = zeros(length(t))

    if !isa(t, Array)
        t = collect(t);
    end

    pos_deb = argmin((t .- td).^2.)

    t_marteau = t[pos_deb:end] .- td

    Ft[pos_deb:end] = F₀*t_marteau.^(k - 1.).*exp.(-t_marteau./θ)/((k - 1.)*θ)^(k - 1.)/exp(1. - k)

    return Ft
end

"""
    smooth_rect(F₀, td, tm, duree, t)

Calcul du vecteur d'excitation pour un signal créneau dont les discontinuités
sont adoucies.
L'excitation est modélisée par une fenêtre de Tuckey

# Paramètres
* F0 : Amplitude de l'excitation [N]
* td : Instant de déclenchement de l'effort [s]
* tm : Temps de montée entre 0 et F₀ [s]
* duree : Durée de l'excitation [s]
* t : Discrétisation temporelle [s]

# Sortie
* Ft : Vecteur d'évolution de l'excitation en fonction du temps [N]
"""
function smooth_rect(F₀, td, tm, duree, t)
    Ft = zeros(length(t))

    pos_deb = argmin((t .- td).^2.)
    pos_fin = argmin((t .- td .- duree).^2.)

    Trect = duree - 2tm
    if isless(Trect, 0.)
        error("Il faut que duree >= 2tm")
    end

    pos_rect_deb = argmin((t .- td .- tm).^2.)
    pos_rect_fin = argmin((t .- td .- duree .+ tm).^2.)

    α = 2tm/duree

    t_monte = t[pos_deb:pos_rect_deb] .- td
    Fmonte = 0.5*F₀*(1. .- cos.(2π*t_monte/α/duree))

    t_desc = t[pos_rect_fin:pos_fin] .- td .- duree
    Fdesc = 0.5*F₀*(1. .- cos.(2π*t_desc/α/duree))

    Ft[pos_deb:pos_rect_deb] = Fmonte
    Ft[pos_rect_deb:pos_rect_fin] .= F₀
    Ft[pos_rect_fin:pos_fin] = Fdesc

    return Ft
end