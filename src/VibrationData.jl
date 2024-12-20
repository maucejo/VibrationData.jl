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

include("models/sdof.jl")
include("models/plate.jl")
include("models/bar_rod.jl")
include("models/beam.jl")
include("models/model.jl")
include("solvers/time_solvers.jl")
include("solvers/FRF_solvers.jl")
include("utils/excitation.jl")
include("utils/noise.jl")

end