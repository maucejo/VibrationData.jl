"""
    modal_model(ωₙ, ηₙ)

Calcul des matrices de masses raideur et amortissement généralisé

# Paramètres
    * p : Structure contenant les données relative à la plaque
    * ωₙ : Vecteur des pulsations propres
    * ξₙ : Amortissement modal

# Sorties
    * Kₙ : Matrice de raideur généralisée
    * Mₙ : Matrice de masse généralisées (matrice identidité, car normalisation à la masse)
    * Cₙ : Matrice de masse généralisée
"""
function modal_model(ωₙ, ξₙ)
    return Diagonal(ωₙ.^2), Diagonal(ones(length(ωₙ))), Diagonal(2ξₙ.*ωₙ)
end