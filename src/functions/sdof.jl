"""
Structure containing the data of a sdof system

# Parameters
* m: Mass [kg]
* ω₀: natural angular frequency [rad/s]
* ξ: Damping ratio
"""
@with_kw struct SDOF
    m :: Float64
    ω₀ ::Float64
    ξ ::Float64
end

"""
    free_response(s::SDOF, t, x₀ = 0., v₀ = 1.)

Compute the free response of a single degree of freedom (SDOF) system.

# Inputs
- s: Structure containing the parameters of the SDOF system
- t: A vector of time points at which to evaluate the response
- x₀: Initial displacement (default is 0.0)
- v₀: Initial velocity (default is 1.0)

# Outputs
- rep: The response of the system at the given time points
- env: The envelope of the response
"""
function free_response(s :: SDOF, t, x₀ = 0., v₀ = 1.)
    (; ω₀, ξ) = s
    nt = length(t)

    if ξ < 1.
        Ω₀ = ω₀*√(1 - ξ^2)
        A = x₀
        B = (v₀ + ξ*ω₀*x₀)/Ω₀

        rep = @. exp(-ξ*ω₀*t)*(A*cos(Ω₀*t) + B*sin(Ω₀*t))
        env = @. exp(-ξ*ω₀*t)*√(A^2 + B^2)

    elseif ξ == 1.
        A = x₀
        B = v₀ + ω₀*x₀

        rep = @. (A + B*t)*exp(-ω₀*t)
        env = zeros(nt)
    else
        β = ω₀*√(ξ^2 - 1)
        A = x₀
        B = (v₀ + ξ*ω₀*x₀)/β

        rep = @. exp(-ξ*ω₀*t)*(A*cosh(β*t) + B*sinh(β*t))
        env = zeros(nt)
    end

    return rep, env
end

"""
    forced_response_harmo(s::SDOF, F, ω, t, x₀ = 0., v₀ = 0.)

Computes the forced response of a single degree of freedom (SDOF) system due to a harmonic excitation

# Inputs
- s: Structure containing the parameters of the SDOF system
- F₀: Amplitude of the excitation [N]
- ω: Frequency of the excitation [rad/s]
- t: A vector of time points at which to evaluate the response
- x₀: Initial displacement (default is 0.)
- v₀: Initial velocity (default is 0.)

# Outputs
- x: Response of the system at the given time points
- xh: Homogeneous solution
- xp: Particular solution
- F: Excitation force
"""
function forced_response_harmo(s:: SDOF, F₀, ω, t, x₀ = 0., v₀ = 0.)
    (; m, ω₀, ξ) = s

    A₀ = F₀/m/√((ω₀^2 - ω^2)^2 + (2ξ*ω*ω₀)^2)
    ϕ = atan(2ξ*ω*ω₀, ω₀^2 - ω^2)

    A = x₀ - A₀*cos(ϕ)
    if ξ < 1.
        Ω₀ = ω₀*√(1 - ξ^2)
        B = (v₀ + ξ*ω₀*A - A₀*sin(ϕ))/Ω₀
        xh = @. (A*cos(Ω₀*t) + B*sin(Ω₀*t))*exp(-ξ*ω₀*t)
    elseif B == 1.
        B = v₀ + ω₀*A - A₀*sin(ϕ)
        xh = @. (A + B*t)*exp(-ω₀*t)
    else
        β = ω₀*√(ξ^2 - 1)
        B = (v₀ + ξ*ω₀*A - A₀*sin(ϕ))/β
        xh = @. (A*cosh(β*t) + B*sinh(β*t))*exp(-ξ*ω₀*t)
    end

    xp = A₀*cos.(ω*t .- ϕ)
    x = xh .+ xp

    return x, xh, xp, F*cos.(ω*t)
end

"""
    forced_response_base_harmo(s::SDOF, F, ω, t, x₀ = 0., v₀ = 0.)

Computes the forced response of a single degree of freedom (SDOF) system due to a harmonic base motion

# Inputs
- s: Structure containing the parameters of the SDOF system
- X₀: Amplitude of the base motion [m]
- ω: Frequency of the excitation [rad/s]
- t: A vector of time points at which to evaluate the response
- x₀: Initial displacement (default is 0.)
- v₀: Initial velocity (default is 0.)

# Outputs
- x: Response of the system at the given time points
- xh: Homogeneous solution
- xp: Particular solution
- X: Excitation base motion
"""
function forced_response_base_harmo(s:: SDOF, X₀, ω, t, x₀ = 0., v₀ = 0.)
    (; ω₀, ξ) = s

    A₀ = X₀*√(ω₀^4 + (2ξ*ω*ω₀)^2)/√((ω₀^2 - ω^2)^2 + (2ξ*ω*ω₀)^2)
    ϕ = atan(2ξ*ω*ω₀, ω₀^2 - ω^2) - atan(2ξ*ω, ω₀)

    A = x₀ - A₀*cos(ϕ)
    if ξ < 1.
        Ω₀ = ω₀*√(1 - ξ^2)
        B = (v₀ + ξ*ω₀*A - A₀*sin(ϕ))/Ω₀
        xh = @. (A*cos(Ω₀*t) + B*sin(Ω₀*t))*exp(-ξ*ω₀*t)
    elseif B == 1.
        B = v₀ + ω₀*A - A₀*sin(ϕ)
        xh = @. (A + B*t)*exp(-ω₀*t)
    else
        β = ω₀*√(ξ^2 - 1)
        B = (v₀ + ξ*ω₀*A - A₀*sin(ϕ))/β
        xh = @. (A*cosh(β*t) + B*sinh(β*t))*exp(-ξ*ω₀*t)
    end

    xp = A₀*cos.(ω*t .- ϕ)
    x = xh .+ xp

    return x, xh, xp, X₀*cos.(ω*t)
end

function forced_response_any(s:: SDOF, F, t, x₀ = 0., v₀ = 0.)
    (; m, ω₀, ξ) = s
    # Time step
    Δt = t[2] - t[1]

    # Impulse response
    if ξ < 1.
        Ω₀ = ω₀*√(1 - ξ^2)
        h = @. exp(-ξ*ω₀*t)*sin(Ω₀*t)/m/Ω₀
    elseif ξ == 1.
        h = @. t*exp(-ω₀*t)/m
    else
        β = ω₀*√(ξ^2 - 1)
        h = @. exp(-ξ*ω₀*t)*sinh(β*t)/m/β
    end

    # Free response
    xh = free_response(s, t, x₀, v₀)[1]

    # Duhamel's integral
    xp = Δt*conv(F, h)[1:length(F)]

    return xh .+ xp
end
