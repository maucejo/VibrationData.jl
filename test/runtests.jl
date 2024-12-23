using VibrationData
using TestItems
using TestItemRunner

@run_package_tests

@testitem "Plate modes" begin
    plaq = Plate(0.6, 0.4, 5e-3, 2.1e11, 7800., 0.3)

    ## Définition des fréquences jusqu'à fmax
    fmax = 1e3
    ωₙ = eigval(plaq, fmax)[1]

    @test round(ωₙ[1]/2π, digits = 2) == 111.33
    @test round(ωₙ[end]/2π, digits = 2) == 933.48
end

@testitem "Longitudinal bar modes" begin
    L = 1.
    S = 3e-4
    E = 2.1e11
    ρ = 7800.

    bar = Bar(L, S, E, ρ)

    fmax = 10e3
    ωₙ = eigval(bar, fmax)[1]

    @test round(ωₙ[1]/2π, digits = 2) == 2594.37
    @test round(ωₙ[end]/2π, digits = 2) == 7783.12
end

@testitem "Torsional bar modes" begin
    L = 1.
    I = π*(1e-2^4)/32.
    J = I
    G = 2.1e11/2(1 + 0.33)
    ρ = 7800.

    rod = Rod(L, I, J, G, ρ)

    fmax = 10e3
    ωₙ = eigval(rod, fmax)[1]

    @test round(ωₙ[1]/2π, digits = 2) == 1590.71
    @test round(ωₙ[end]/2π, digits = 2) == 9544.27
end

@testitem "Beam modes" begin
    L = 1.
    b = 3e-2
    h = 1e-2
    S = b*h
    I = b*h^3/12.
    E = 2.1e11
    ρ = 7800.

    beam = Beam(L, S, I, E, ρ)

    fmax = 1e3
    ωₙ = eigval(beam, fmax)[1]

    @test round(ωₙ[1]/2π, digits = 2) == 23.53
    @test round(ωₙ[end]/2π, digits = 2) == 847.02
end

@testitem "Rectangular shape force" begin
    Δt = 1e-6 # Pas de temps
    tf = 0.07 # instant final
    t = 0.:Δt:tf

    rect = Rectangle(1., 8e-3, 1e-2)
    Ft = excitation(rect, t)

    pos = findall(Ft .== 1.)
    @test maximum(Ft) == 1.
    @test sum(diff(Ft)) == 0.
end

@testitem "Impact force" begin
    Δt = 1e-6 # Pas de temps
    tf = 0.07 # instant final
    t = 0.:Δt:tf

    hammer = Hammer(1., 8e-3, 9.7, 6e-4)
    Ft = excitation(hammer, t)

    @test round(maximum(Ft), digits = 2) == 1.
    @test round(length(t)*sum(Ft)*Δt) == 314.
end

@testitem "Triangle shape force" begin
    Δt = 1e-6 # Pas de temps
    tf = 0.07 # instant final
    t = 0.:Δt:tf

    triangle = Triangle(1. , 8e-3, 5e-2)
    Ft = excitation(triangle, t)

    @test round(maximum(Ft), digits = 2) == 1.
    @test isapprox(sum(diff(Ft)), 0., atol = eps())
end

@testitem "Smooth rectangle shape force" begin
    Δt = 1e-6 # Pas de temps
    tf = 0.07 # instant final
    t = 0.:Δt:tf

    srect = SmoothRect(1., 8e-3, 5e-3, 5e-2)
    Ft = excitation(srect, t)

    @test round(maximum(Ft), digits = 2) == 1.
    @test sum(diff(Ft)) == 0.
end

@testitem "Random force" begin
    Δt = 1e-6 # Pas de temps
    tf = 0.07 # instant final
    t = 0.:Δt:tf

    randexc = RandomExc(1., 8e-3, 5e-2, 0.1)
    Ft = excitation(randexc, t)

    Frand = Ft[8e-3 .≤ t .≤ 5.8e-2]
    Fmean = sum(Frand)/length(Frand)
    Fstd = sqrt(sum(abs2, Frand .- Fmean)/(length(Frand) - 1))
    @test abs(round(Fmean, digits = 2) - 1.) ≤ 1e-2
    @test abs(round(Fstd, digits = 2) - 0.1)/0.1 ≤ 1e-2
end

@testitem "Time solvers" begin
    plate = Plate(0.6, 0.4, 5e-3, 2.1e11, 7800., 0.3)

    fmax = 1e3
    ωₙ, kₙ = eigval(plate, fmax)
    Nmodes = length(ωₙ)

    Kₙ, Mₙ, Cₙ = modal_model(ωₙ, 1e-2)

    Δt = 1e-6 # Pas de temps
    tf = 0.07 # instant final
    t = 0.:Δt:tf
    loc = [0.1, 0.2]

    hammer = Hammer(1., 8e-3, 9.7, 6e-4)
    F = excitation(hammer, t)
    ϕₑ = eigmode(plate, kₙ, loc[1], loc[2])
    Fₙ = (F*ϕₑ)'

    ϕₒ = eigmode(plate, kₙ, loc[1], loc[2])

    prob = LinearTimeProblem(Kₙ, Mₙ, Cₙ, Fₙ, t)
    CI = (D₀ = zeros(Nmodes), V₀ = zeros(Nmodes))

    # Generalized-α
    solGα = solve(prob, CI, GeneralizedAlpha())
    (; A) = solGα
    AccGα = ϕₒ*A

    # Central difference
    solCD = solve(prob, CI, CentralDiff())
    (; A) = solCD
    AccCD = ϕₒ*A

    # HHT
    solHHT = solve(prob, CI, HHT())
    (; A) = solHHT
    AccHHT = ϕₒ*A

    # Fox-Goodwin
    solFG = solve(prob, CI, FoxGoodwin())
    (; A) = solFG
    AccFG = ϕₒ*A

    # Linear acceleration
    solLA = solve(prob, CI, LinearAcceleration())
    (; A) = solLA
    AccLA = ϕₒ*A

    # Newmark
    solNM = solve(prob, CI, Newmark())
    (; A) = solNM
    AccNM = ϕₒ*A

    # WBZ
    solWBZ = solve(prob, CI, WBZ())
    (; A) = solWBZ
    AccWBZ = ϕₒ*A

    # Mid-point
    solMP = solve(prob, CI, MidPoint())
    (; A) = solMP
    AccMP = ϕₒ*A

    # RK4
    solRK4 = solve(prob, CI, RK4())
    (; A) = solRK4
    AccRK4 = ϕₒ*A

    energyGα = sum(abs2, AccGα)
    energyCD = sum(abs2, AccCD)
    energyHHT = sum(abs2, AccHHT)
    energyFG = sum(abs2, AccFG)
    energyLA = sum(abs2, AccLA)
    energyNM = sum(abs2, AccNM)
    energyWBZ = sum(abs2, AccWBZ)
    energyMP = sum(abs2, AccMP)
    energyRK4 = sum(abs2, AccRK4)

    @test round(sum(abs2, AccGα), digits = 2) == 866.97
    @test (abs(energyCD - energyGα)/energyGα) ≤ 1e-2
    @test (abs(energyHHT - energyGα)/energyGα) ≤ 1e-2
    @test (abs(energyFG - energyGα)/energyGα) ≤ 1e-2
    @test (abs(energyLA - energyGα)/energyGα) ≤ 1e-2
    @test (abs(energyNM - energyGα)/energyGα) ≤ 1e-2
    @test (abs(energyWBZ - energyGα)/energyGα) ≤ 1e-2
    @test (abs(energyMP - energyGα)/energyGα) ≤ 1e-2
    @test (abs(energyRK4 - energyGα)/energyGα) ≤ 1e-2
end

@testitem "Modal FRF" begin
    plaq = Plate(0.6, 0.4, 5e-3, 2.1e11, 7800., 0.3)

    ## Définition des fréquences jusqu'à fmax
    fmax = 1e3
    ωₙ, kₙ = eigval(plaq, fmax)
    Nmodes = length(ωₙ)

    ## Définition du modèle modal
    Kₙ, Mₙ, Cₙ = modal_model(ωₙ, 1e-2)

    # Définition de la déformée modale
    loc = [0.1, 0.1]
    ϕₑ = eigmode(plaq, kₙ, loc[1], loc[2])

    # Calcul des coordonnées généralisées
    freq = 70:125
    ξ = 1e-2
    prob = ModalFRF(ωₙ, ξ, ϕₑ, ϕₑ, freq)

    _FRF = frf(prob, :acc)

    FRF = reduce(hcat, _FRF)[:]

    f = ωₙ[1]/2π
    pos_max = argmax(abs.(FRF))
    @test (abs(freq[pos_max] - f))/f ≤ 1e-2

    maxAccth = ϕₑ[:, 1][1]^2/2ξ
    maxAcc = maximum(abs.(FRF))
    @test abs(maxAcc - maxAccth)/maxAccth ≤ 3e-2
end