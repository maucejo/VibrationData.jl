"""
    excitation(param, t)

Computes different types of excitation signals

# Parameters
* param : NamedTuple or struct containing the following fields:
    * type : type of excitation
        1. :triangle
        2. :rectangle
        3. :marteau
        4. :srect (smooth rectangle)
    * F₀ : Amplitude of the force [N]
    * td : Time of force initiation [s]
    * duree : Duration of the excitation [s]
        *note : Mandatory for types 1, 2, and 4, not needed for type 3*
    * k and θ : Shape parameter [-] and intensity parameter [s]
        *note : Mandatory for type 3, not needed for other types*
    * tm : Rise time from 0 to F₀ [s]
        *note : Mandatory for type 4, not needed for other types*

# Output
* F : Vector of excitation evolution over time [N]
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

Computes the excitation vector for a triangular signal

# Parameters
* F₀ : Amplitude of the excitation [N]
* td : Time of force initiation [s]
* duree : Duration of the excitation [s]
* t : Time discretization [s]

# Output
* Ft : Vector of excitation evolution over time [N]
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

Computes the excitation vector for a rectangular signal

# Parameters
* F₀ : Amplitude of the excitation [N]
* td : Time of force initiation [s]
* duree : Duration of the excitation [s]
* t : Time discretization [s]

# Output
* Ft : Vector of excitation evolution over time [N]
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

Computes the excitation vector for a hammer impact signal
The hammer impact is assumed to have the shape of a gamma distribution
with parameters (k, theta)

# Parameters
* F₀ : Amplitude of the excitation [N]
* td : Time of force initiation [s]
* k : Shape parameter
* θ : Intensity parameter [s]
* t : Time discretization [s]

# Output
* Ft : Vector of excitation evolution over time [N]
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

Computes the excitation vector for a smooth rectangular signal.
The excitation is modeled by a Tuckey window.

# Parameters
* F₀ : Amplitude of the excitation [N]
* td : Time of force initiation [s]
* tm : Rise time from 0 to F₀ [s]
* duree : Duration of the excitation [s]
* t : Time discretization [s]

# Output
* Ft : Vector of excitation evolution over time [N]
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