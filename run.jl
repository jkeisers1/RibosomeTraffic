using Pkg
Pkg.activate(@__DIR__)

using RibosomeTraffic
using Printf
using BenchmarkTools
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

k_pause = 0.001      # Small pause rate to test the stats
k_unpause = 0.1      # Recovery rate
lifetime_mRNA = 300.0
delta = 1.0 / lifetime_mRNA

model = TranscriptModel(L, l_ribosome, α, β, k_elong, k_pause, k_unpause, delta)

# ==========================================
# 2. SIMULATION
# ==========================================
t_bench = 100.0
println("--- Running Simulation (t_max = $t_bench s) ---")

@btime state = run_custom_simulation(model, t_bench)

# ==========================================
# 3. ANALYSIS
# ==========================================

# --- A. Particle Types (Microscopic) ---
# "What are the particles doing?"
avg_N_active = state.cum_active_time / state.time
avg_N_paused = state.cum_paused_time / state.time
avg_N_mobile = state.cum_mobile_time / state.time
avg_N_jammed = avg_N_active - avg_N_mobile

rho_total  = (avg_N_active + avg_N_paused) / L
rho_active = avg_N_active / L
rho_paused = avg_N_paused / L
rho_mobile = avg_N_mobile / L
rho_jammed = avg_N_jammed / L

J_total = state.flux_termination / state.time

# --- B. System States (Macroscopic) ---
# "Is the whole mRNA clean or dirty?"

# 1. Time Split
time_unpaused = state.cum_time_unpaused
time_paused   = state.time - state.cum_time_unpaused
frac_unpaused = time_unpaused / state.time

# 2. Flux Split (Avoid division by zero)
J_unpaused = (time_unpaused > 0) ? (state.flux_unpaused / time_unpaused) : 0.0
J_paused   = (time_paused > 0)   ? (state.flux_paused / time_paused) : 0.0

# 3. Density Split (Total density in that state)
rho_sys_unpaused = (time_unpaused > 0) ? (state.cum_mass_unpaused / time_unpaused / L) : 0.0
rho_sys_paused   = (time_paused > 0)   ? (state.cum_mass_paused / time_paused / L) : 0.0

# ==========================================
# 4. REPORT
# ==========================================
println("\n================ RESULTS ================")
@printf "Total Time:       %.2f s\n" state.time
@printf "Total Flux (J):   %.4f ribosomes/sec\n" J_total
@printf "Total Density:    %.4f /codon\n" rho_total
println("-----------------------------------------")
println("PARTICLE BREAKDOWN (Microscopic)")
println("  ├─ Active Density:  $(round(rho_active, digits=4))")
println("  │   ├─ Mobile:      $(round(rho_mobile, digits=4))   (Free to move)")
println("  │   └─ Jammed:      $(round(rho_jammed, digits=4))   (Blocked by neighbor)")
println("  └─ Paused Density:  $(round(rho_paused, digits=4))   (Chemically paused)")
println("-----------------------------------------")
println("SYSTEM STATE BREAKDOWN (Macroscopic)")
@printf "  [Unpaused State] (%.1f%% of time)\n" (frac_unpaused * 100)
@printf "   ├─ Flux:    %.4f\n" J_unpaused
@printf "   └─ Density: %.4f\n" rho_sys_unpaused
println()
@printf "  [Paused State]   (%.1f%% of time)\n" ((1.0 - frac_unpaused) * 100)
@printf "   ├─ Flux:    %.4f\n" J_paused
@printf "   └─ Density: %.4f\n" rho_sys_paused
println("=========================================")