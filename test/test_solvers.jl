using LazyGrids
using VibrationData

## Définition de la structure - Plaque
plaq = Plate(0.6, 0.4, 5e-3, 2.1e11, 7800., 0.3)

## Définition des fréquences jusqu'à fmax
fmax = 1e3
ωₙ, kₙ = eigval(plaq, 10fmax)
Nmodes = length(ωₙ)

## Définition du modèle modal
Kₙ, Mₙ, Cₙ = modal_model(ωₙ, 1e-2)

#%% Définition de l'effort généralisé
Δt = 1e-6 # Pas de temps
tf = 0.07 # instant final
t = 0.:Δt:tf
loc = [0.1, 0.2]

param = (type = :marteau, F₀ = 1., td = 8e-3, k = 9.7, θ = 6e-4)

ϕₑ = eigmode(plaq, kₙ, loc[1], loc[2])
F = excitation(param, t)
Fₙ = (F*ϕₑ)'

## Définition des déformées modales aux points d'écoute
Nx = 40
Ny = 40
x, y = ndgrid(LinRange(0., plaq.L, Nx), LinRange(0., plaq.b, Ny))
X = x[:]
Y = y[:]
ϕₒ = eigmode(plaq, kₙ, X, Y)

# Calcul des coordonnées généralisées
prob = LinearTimeProblem(Kₙ, Mₙ, Cₙ, Fₙ, t)
u0 = (D₀ = zeros(Nmodes), V₀ = zeros(Nmodes))

sol = solve(prob, u0)
(; D, V, A) = sol

Dep = ϕₒ*D
Vit = ϕₒ*V
Acc = ϕₒ*A