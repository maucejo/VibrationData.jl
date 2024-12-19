"""
    excitation(param, t)

Computes different types of excitation signals

# Parameters
* param : NamedTuple or struct containing the following fields:
    * type : type of excitation
        1. :triangle
        2. :rectangle
        3. :hammer
        4. :srect (smooth rectangle)
    * F₀ : Amplitude of the force [N]
    * ts : Time of force initiation [s]
    * duration : Duration of the excitation [s]
        *note : Mandatory for types 1, 2, and 4, not needed for type 3*
    * k and θ : Shape parameter [-] and intensity parameter [s]
        *note : Mandatory for type 3, not needed for other types*
    * tr : Rise time from 0 to F₀ [s]
        *note : Mandatory for type 4, not needed for other types*

# Output
* F : Vector of excitation evolution over time [N]
"""
function excitation(param, t)
    (; type, F₀, ts) = param

    if type == :triangle
        (; duration) = param
        F = triangle(F₀, ts, duration, t)
    elseif type == :rectangle
        (; duration) = param
        F = rectangle(F₀, ts, duration, t)
    elseif type == :random
        (; duration, σ)
        F = random(F₀, ts, duration, σ, t)
    elseif type == :hammer
        (; k, θ) = param
        F = hammer(F₀, ts, k, θ, t)
    elseif type == :srect
        (; tr, duration) = param
        F = smooth_rect(F₀, ts, tr, duration, t)
    else
        error("Excitation type not implemented")
    end
end

"""
    triangle(F₀, ts, duration, t)

Computes the excitation vector for a triangular signal

# Parameters
* F₀ : Amplitude of the excitation [N]
* ts : Time of force initiation [s]
* duration : Duration of the excitation [s]
* t : Time discretization [s]

# Output
* Ft : Vector of excitation evolution over time [N]
"""
function triangle(F₀, ts, duration, t)
    Ft = zeros(length(t))

    tr = (2ts + duration)/2.
    pos_start = argmin((t .- ts).^2.)
    pos_middle = argmin((t .- tr).^2.)
    pos_end = argmin((t .- ts .- duration).^2.)
    amp = 2F₀/duration

    Ft[pos_start:pos_middle] = amp*(t[pos_start:pos_middle] .- ts)
    Ft[pos_middle + 1:pos_end] = F₀ .- amp*(t[pos_middle + 1:pos_end] .- tr)

    return Ft
end

"""
    rectangle(F₀, ts, duration, t)

Computes the excitation vector for a rectangular signal

# Parameters
* F₀ : Amplitude of the excitation [N]
* ts : Time of force initiation [s]
* duration : Duration of the excitation [s]
* t : Time discretization [s]

# Output
* Ft : Vector of excitation evolution over time [N]
"""
function rectangle(F₀, ts, duration, t)
    Ft = zeros(length(t))

    pos_start = argmin((t .- ts).^2.)
    pos_end = argmin((t .- ts .- duration).^2.)

    pos_exc_t = findall(t[pos_start] .≤ t .≤ t[pos_end])

    Ft[pos_exc_t] .= F₀

    return Ft
end

"""
    random(F₀, ts, duration, σ, t)

Computes the excitation vector for a random signal

# Parameters
* F₀ : Amplitude of the excitation [N]
* ts : Time of force initiation [s]
* duration : Duration of the excitation [s]
* σ : Standard deviation of the random noise
* t : Time discretization [s]

# Output
* Ft : Vector of excitation evolution over time [N]
"""
function random(F₀, ts, duration, σ, t)
    Ft = zeros(length(t))

    pos_start = argmin((t .- ts).^2.)
    pos_end = argmin((t .- ts .- duration).^2.)

    pos_exc_t = findall(t[pos_start] .≤ t .≤ t[pos_end])

    Ft[pos_exc_t] .= F₀ + σ*randn(length(pos_exc_t))

    return Ft
end

"""
    hammer(F₀, ts, k, θ, t)

Computes the excitation vector for a hammer impact signal
The hammer impact is assumed to have the shape of a gamma distribution
with parameters (k, theta)

# Parameters
* F₀ : Amplitude of the excitation [N]
* ts : Time of force initiation [s]
* k : Shape parameter
* θ : Intensity parameter [s]
* t : Time discretization [s]

# Output
* Ft : Vector of excitation evolution over time [N]
"""
function hammer(F₀, ts, k, θ, t)
    Ft = zeros(length(t))

    if !isa(t, Array)
        t = collect(t);
    end

    pos_start = argmin((t .- ts).^2.)

    t_hammer = t[pos_start:end] .- ts

    Ft[pos_start:end] = F₀*t_hammer.^(k - 1.).*exp.(-t_hammer./θ)/((k - 1.)*θ)^(k - 1.)/exp(1. - k)

    return Ft
end

"""
    smooth_rect(F₀, ts, tr, duration, t)

Computes the excitation vector for a smooth rectangular signal.
The excitation is modeled by a Tuckey window.

# Parameters
* F₀ : Amplitude of the excitation [N]
* ts : Time of force initiation [s]
* tr : Rise time from 0 to F₀ [s]
* duration : Duration of the excitation [s]
* t : Time discretization [s]

# Output
* Ft : Vector of excitation evolution over time [N]
"""
function smooth_rect(F₀, ts, tr, duration, t)
    Ft = zeros(length(t))

    pos_start = argmin((t .- ts).^2.)
    pos_end = argmin((t .- ts .- duration).^2.)

    Trect = duration - 2tr
    if isless(Trect, 0.)
        error("Il faut que duration >= 2tr")
    end

    pos_rect_start = argmin((t .- ts .- tr).^2.)
    pos_rect_end = argmin((t .- ts .- duration .+ tr).^2.)

    α = 2tr/duration

    t_rise = t[pos_start:pos_rect_start] .- ts
    Frise = 0.5*F₀*(1. .- cos.(2π*t_rise/α/duration))

    t_desc = t[pos_rect_end:pos_end] .- ts .- duration
    Fdesc = 0.5*F₀*(1. .- cos.(2π*t_desc/α/duration))

    Ft[pos_start:pos_rect_start] = Frise
    Ft[pos_rect_start:pos_rect_end] .= F₀
    Ft[pos_rect_end:pos_end] = Fdesc

    return Ft
end