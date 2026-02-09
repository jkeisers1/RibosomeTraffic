using RibosomeTraffic
using Test
using Statistics
using Pkg
# If we are running this as a script, ensure the current project is active
if !("RibosomeTraffic" in keys(Pkg.project().dependencies))
    Pkg.activate(joinpath(@__DIR__, ".."))
end


println("--- Starting Regulatory-Grade Unit Tests ---")


@testset "Solver Equivalence Tests" begin
    # 1. Define a Small, Fast Model for Testing
    # We use a small L so the SciML solver doesn't take forever.
    L = 50
    l_rib = 5
    model = TranscriptModel(
        L,
        l_rib,
        0.3,    # alpha (initiation)
        0.5,    # beta (termination)
        ones(L),# k_elong
        0.01,   # k_pause
        0.1,    # k_unpause
        0.0,     # no degradation for this test
    )

    t_max = 500.0
    N_trials = 50 # Number of repeats to smooth out noise

    println("   Running comparison on L=$L lattice ($N_trials trials)...")

    # --- 2. Run the Golden Master (SciML) ---
    # We collect the final number of ribosomes for each run
    sciml_results = zeros(N_trials)
    for i = 1:N_trials
        # We use the helper we wrote in solver_jump.jl
        sol = run_sciml_simulation(model, t_max)

        # Calculate average density from the solution
        # sol.u is a vector of arrays.
        # We just take the final state for simplicity here.
        final_lattice = sol.u[end]
        sciml_results[i] = count(x -> x > 0, final_lattice)
    end
    mean_sciml = mean(sciml_results)
    println("   > Mean Ribosomes (SciML):  $mean_sciml")

    # --- 3. Run the Candidate (Custom) ---
    custom_results = zeros(N_trials)
    for i = 1:N_trials
        state = run_custom_simulation(model, t_max)

        # Our custom state tracks counts automatically!
        # We take the time-averaged density (even better than final state)
        # avg_N = (cum_active + cum_paused) / time
        avg_particles = (state.cum_active_time + state.cum_paused_time) / state.time
        custom_results[i] = avg_particles
    end
    mean_custom = mean(custom_results)
    println("   > Mean Ribosomes (Custom): $mean_custom")

    # --- 4. The "Regulatory" Check ---
    # Are they within 10% of each other?
    # (Stochastic noise means they won't be identical, but they should be close)
    diff_percent = abs(mean_sciml - mean_custom) / mean_sciml * 100
    println("   > Difference: $(round(diff_percent, digits=2))%")

    @test diff_percent < 10.0
end
