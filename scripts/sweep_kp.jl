# scripts/3_sweep_kp.jl
using Distributed
using CSV, DataFrames  # You'll need: add CSV DataFrames
using Statistics
# 1. Setup Parallelism
if nprocs() == 1
    addprocs(4) # Adjust for your machine
end

@everywhere begin
    using Pkg; Pkg.activate(".")
    using RibosomeTraffic
    
    # Define the worker function inside here so workers can see it
    function run_point(kp_val)
        # Create model
        L = 1000
        # ... standard params ...
        model = TranscriptModel(L, 10, 0.1, 2.0, ones(L), kp_val, 0.1, 0.0033)
        
        # Run Long Simulation
        state = run_custom_simulation(model, 50000.0)
        
        # Return a NamedTuple (easy to turn into a DataFrame later)
        return (
            kp = kp_val,
            flux  = state.flux_termination / state.time,
            density = (state.cum_active_time / state.time) / L,
            jammed  = ((state.cum_active_time - state.cum_moving_masses)/state.time) / L
        )
    end
end

# 2. Define the List (The Sweep)
# Run 10 values, 3 repeats each = 30 simulations
kp_list = range(0.1, 1.2, length=10)
work_list = repeat(kp_list, inner=3)

println("--- Starting Sweep of $(length(work_list)) simulations ---")

# 3. The Parallel Map
results = pmap(run_point, work_list)

# 4. Save to Disk
# Convert the list of NamedTuples into a Table
df = DataFrame(results)

# Calculate mean per kp (Optional, or do this in your plotting script)
gdf = combine(groupby(df, :kp), :flux => mean, :density => mean)

println("--- Saving Results ---")
CSV.write("results_sweep_kp.csv", df)
CSV.write("results_summary.csv", gdf)
println("Done! Data saved to CSV.")