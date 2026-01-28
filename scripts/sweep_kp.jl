using Distributed
using CSV, DataFrames
using Statistics
using UUIDs
using Dates
using Pkg # Needed for environment setup

# ==============================================================================
#  SECTION 0: ENVIRONMENT SETUP (Crucial Step!)
# ==============================================================================

# 1. Define the absolute path to the project root (one level up from scripts/)
project_root = abspath(joinpath(@__DIR__, ".."))

println("--- Activating Environment at: $project_root ---")

# 2. Activate for the MAIN process (The Manager)
Pkg.activate(project_root)
using RibosomeTraffic # Now Main knows what TranscriptModel is

# ==============================================================================
#  SECTION 1: EXPERIMENT CONFIGURATION
# ==============================================================================
CONSTANTS = Dict(
    :L          => 1000,
    :l_rib      => 10,
    :alpha      => 0.6,
    :beta       => 1.0,
    :k_unpause  => 0.1,
    :lifetime   => 300.0,
    :t_max      => 10^5,
    :repeats    => 4
)

sweep_variable = "k_pause"
sweep_values   = collect(range(0.0, 0.5, length=12))

# ==============================================================================
#  SECTION 2: AUTOMATIC LOGGING GENERATOR
# ==============================================================================
batch_id = string(uuid4())[1:8]
timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM")
desc = "Sweep_$(sweep_variable)_$(first(sweep_values))_to_$(last(sweep_values))_a$(CONSTANTS[:alpha])"

println("\n╔══════════════════════════════════════════════════════════╗")
println("║  🧪 STARTING EXPERIMENT: $batch_id                ")
println("║  📅 Date: $timestamp                                     ")
println("║  📝 Description: $desc                                   ")
println("╚══════════════════════════════════════════════════════════╝")

log_entry = DataFrame(
    batch_id      = batch_id,
    date          = timestamp,
    description   = desc,
    sweep_var     = sweep_variable,
    sweep_range   = string(first(sweep_values)) * " -> " * string(last(sweep_values)),
    L             = CONSTANTS[:L],
    alpha         = CONSTANTS[:alpha],
    beta          = CONSTANTS[:beta],
    k_unpause     = CONSTANTS[:k_unpause],
    lifetime      = CONSTANTS[:lifetime],
    repeats       = CONSTANTS[:repeats]
)

# ==============================================================================
#  SECTION 3: PARALLEL EXECUTION
# ==============================================================================
if nprocs() == 1
    addprocs(max(1, Sys.CPU_THREADS - 2)) 
    println("--- Added $(nworkers()) worker processes ---")
end

# 3. Activate for WORKER processes
# We use $project_root to pass the EXACT string calculated above
@everywhere begin
    using Pkg
    Pkg.activate($project_root) 
    using RibosomeTraffic
    
    function run_point(args)
        val, params = args 
        L = params[:L]
        delta = 1.0 / params[:lifetime]
        
        # --- MODEL SETUP ---
        current_kp = val 
        model = TranscriptModel(L, params[:l_rib], params[:alpha], params[:beta], 
                                ones(L), current_kp, params[:k_unpause], delta)
        
        # --- RUN SIMULATION ---
        state = run_custom_simulation(model, params[:t_max])
        
        # --- CALCULATE OBSERVABLES ---
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

        time_unpaused = state.cum_time_unpaused
        time_paused   = state.time - state.cum_time_unpaused
        frac_unpaused = time_unpaused / state.time

        J_unpaused = (time_unpaused > 0) ? (state.flux_unpaused / time_unpaused) : 0.0
        J_paused   = (time_paused > 0)   ? (state.flux_paused / time_paused) : 0.0

        rho_sys_unpaused = (time_unpaused > 0) ? (state.cum_mass_unpaused / time_unpaused / L) : 0.0
        rho_sys_paused   = (time_paused > 0)   ? (state.cum_mass_paused / time_paused / L) : 0.0
        
        return (
            param_val = val,
            J_total    = J_total,
            rho_total  = rho_total,
            rho_active = rho_active,
            rho_paused = rho_paused,
            rho_mobile = rho_mobile,
            rho_jammed = rho_jammed,
            frac_unpaused    = frac_unpaused,
            J_unpaused       = J_unpaused,
            J_paused         = J_paused,
            rho_sys_unpaused = rho_sys_unpaused,
            rho_sys_paused   = rho_sys_paused
        )
    end
end

work_list = []
for v in sweep_values
    for r in 1:CONSTANTS[:repeats]
        push!(work_list, (v, CONSTANTS))
    end
end

println("--- Launching $(length(work_list)) simulations... ---")
@time results = pmap(run_point, work_list)

# ==============================================================================
#  SECTION 4: DATA PROCESSING & SAVING
# ==============================================================================
df_raw = DataFrame(results)
df_raw[!, :batch_id] .= batch_id 

# Summary Stats
gdf = combine(groupby(df_raw, :param_val), 
    Cols(Not([:param_val, :batch_id])) .=> mean
)
gdf[!, :batch_id] .= batch_id

if !isdir("data"); mkdir("data"); end
log_file = "experiment_log.csv"

println("--- Saving Data ---")

if !isfile(log_file)
    CSV.write(log_file, log_entry)
else
    CSV.write(log_file, log_entry, append=true)
end

CSV.write("data/$(batch_id)_raw.csv", df_raw)
CSV.write("data/$(batch_id)_summary.csv", gdf)

println("✅ SUCCESS! Batch ID: $batch_id")
println("   📁 Log:        $log_file")
println("   📊 Summary:    data/$(batch_id)_summary.csv")