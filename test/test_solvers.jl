using LazyGrids
using VibrationData

## Plate definition
plate = Plate(0.6, 0.4, 5e-3, 2.1e11, 7800., 0.3)

## Computation of the resonance frequencies
fmax = 1e3
ωₙ, kₙ = eigval(plate, 10fmax)
Nmodes = length(ωₙ)

## Construction of the modale model
Kₙ, Mₙ, Cₙ = modal_model(ωₙ, 1e-2)

#%% Calculation of the modal force vector
Δt = 1e-6 # Pas de temps
tf = 0.07 # instant final
t = 0.:Δt:tf
loc = [0.1, 0.2]

hammer = Hammer(1., 8e-3, 9.7, 6e-4)
F = excitation(hammer, t)
ϕₑ = eigmode(plate, kₙ, loc[1], loc[2])
Fₙ = (F*ϕₑ)'

## Listening points mesh
Nx = 40
Ny = 40
x, y = ndgrid(LinRange(0., plate.L, Nx), LinRange(0., plate.b, Ny))
X = x[:]
Y = y[:]
ϕₒ = eigmode(plate, kₙ, X, Y)

# Computation of the modal coordinates
prob = LinearTimeProblem(Kₙ, Mₙ, Cₙ, Fₙ, t)
u0 = (D₀ = zeros(Nmodes), V₀ = zeros(Nmodes))

sol = solve(prob, u0)
(; D, V, A) = sol

## Computation of the displacement, velocity and acceleration at the listening points
Disp = ϕₒ*D
Vel = ϕₒ*V
Acc = ϕₒ*A