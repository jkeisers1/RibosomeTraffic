module RibosomeTraffic

# --- Dependencies ---
using StaticArrays
using Distributions
using JumpProcesses
using Random

# --- Exported Interface ---
# These are the functions users (and you) will call:
export TranscriptModel, SimState      # From types.jl
export step!, run_custom_simulation   # From solver_custom.jl
export build_jump_problem             # From solver_jump.jl
export plot_kymograph                 # From utils.jl

# --- Include Sub-Modules ---
# The order matters! Types must be loaded first.

# 1. Define the Data Structures
include("types.jl")

# 2. Load the Solvers
include("solver_custom.jl") # Your high-perf custom loop
include("solver_jump.jl")   # The SciML integration

# 3. Load Utilities (Analysis/Plotting helpers)
include("utils.jl")

end # module
