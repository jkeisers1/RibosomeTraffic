using Pkg
Pkg.activate(@__DIR__)

using BenchmarkTools
using RibosomeTraffic
using Printf

# ==========================================
# 1. SETUP
# ==========================================
println("--- Initializing Model ---")

L = 1000             # Length of transcript
l_ribosome = 10      # Ribosome footprint
α = 0.6              # Initiation rate
β = 1.0              # Termination rate

# Create a heterogeneous elongation profile
k_elong = ones(L)

k_pause = 0.0     # Pause probability
k_unpause = 0.0     # Recovery rate
tau = 120
delta = 1.0 / tau  # mRNA degradation rate
model = TranscriptModel(L, l_ribosome, α, β, k_elong, k_pause, k_unpause, delta)

# ==========================================
# 2. SIMULATION
# ==========================================
#t_max = 5000.0
t_bench = 1000.0
println("--- Running Simulation (t_max = \$t_max s) ---")

@time state = run_custom_simulation(model, t_bench)

# ==========================================
# 3. ANALYSIS (Time-Averaged)
# ==========================================

# A. Average Particle Counts (N)
# We divide the integrated "Particle-Seconds" by the Total Time
avg_N_active = state.cum_active_time / state.time
avg_N_paused = state.cum_paused_time / state.time
avg_N_moving = state.cum_moving_masses / state.time # Effectively "N_free"

# Derived: Jammed = Active (Total) - Moving (Free)
avg_N_jammed = avg_N_active - avg_N_moving

# B. Average Densities (ρ = N / L)
rho_total  = (avg_N_active + avg_N_paused) / L
rho_active = avg_N_active / L
rho_paused = avg_N_paused / L
rho_jammed = avg_N_jammed / L

# C. Flux
J = state.flux_termination / state.time

# ==========================================
# 4. REPORT
# ==========================================
println("\n--- Results (Time-Averaged) ---")
@printf "Total Time:       %.2f s\n" state.time
@printf "Total Steps:      %d\n" state.step_count
@printf "Flux (J):         %.4f ribosomes/sec\n" J
println("--------------------------------")
@printf "Total Density:    %.4f /codon \n" rho_total 
println("  ├─ Active:      $(round(rho_active, digits=4))")
println("  │   ├─ Moving:  $(round(avg_N_moving / L, digits=4))")
println("  │   └─ Jammed:  $(round(rho_jammed, digits=4))")
println("  └─ Paused:      $(round(rho_paused, digits=4))")
println("--------------------------------")
