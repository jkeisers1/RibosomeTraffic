using Pkg
# Ensure we use the environment in this current folder
Pkg.activate(@__DIR__)

using RibosomeTraffic
using Printf

# ==========================================
# 1. SETUP
# ==========================================
println("--- Initializing Model ---")

# Biological Parameters
L = 1000             # Length of transcript
l_ribosome = 1     # Ribosome footprint
α = 0.2              # Initiation rate (0.6/sec)
β = 1.0              # Termination rate
k_pause = 0.0      # Pause probability
k_unpause = 0.0      # Recovery rate

# Elongation Profile (Uniform for now)
k_elong = ones(L)

# Initialize the Model
model = TranscriptModel(L, l_ribosome, α, β, k_elong, k_pause, k_unpause)

# ==========================================
# 2. SIMULATION
# ==========================================
t_max = 100000.0
println("--- Running Simulation (t_max = $t_max s) ---")

# The @time macro will print: seconds taken, allocations, etc.
@time state = run_custom_simulation(model, t_max)

# ==========================================
# 3. ANALYSIS
# ==========================================
J = state.flux_termination / state.time
rho = count(x -> x == 1, state.lattice) / L

println("\n--- Results ---")
@printf "Total Time:       %.2f s\n" state.time
@printf "Total Steps:      %d\n" state.step_count
@printf "Protein Output:   %d\n" state.flux_termination
@printf "Flux (J):         %.4f ribosomes/sec\n" J
@printf "Density (ρ):      %.4f ribosomes/codon\n" rho
println("--------------------------------")