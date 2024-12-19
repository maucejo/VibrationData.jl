module VibrationData

using Parameters, ProgressMeter, LinearAlgebra, DSP

# Structs
export Plate, Bar, Rod, Beam, SDOF,
       LinearTimeProblem, TimeSolution,
       ModalFRF, DirectFRF

# Functions
export excitation,
       eigval, eigmode, modal_model,
       solve,
       frf

# Time solvers
export CentralDiff, RK4, FoxGoodwin, LinearAcceleration,
       Newmark, HHT, WBZ, GeneralizedAlpha, MidPoint

# Noise utils
export agwn, varest, estimated_SNR

include("functions/excitation.jl")
include("functions/sdof.jl")
include("functions/plate.jl")
include("functions/bar_rod.jl")
include("functions/beam.jl")
include("functions/model.jl")
include("functions/time_solvers.jl")
include("functions/FRF_solvers.jl")
include("functions/noise.jl")

end