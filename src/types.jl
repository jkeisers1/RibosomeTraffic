struct TranscriptModel
    L::Int
    l_ribosome::Int
    α::Float64              # Initiation rate
    β::Float64              # Termination rate
    k_elong::Vector{Float64}
    k_pause::Float64
    k_unpause::Float64
    delta::Float64          # Degradation rate
end

mutable struct SimState
    time::Float64
    step_count::Int
    
    # Lattice & State
    lattice::Vector{Int} 
    internal_states::Vector{UInt8} # 1=Active, 2=Paused
    
    # Propensities
    rate_elong::Vector{Float64}
    rate_initiation::Float64
    rate_switch::Vector{Float64}
    total_rate_elong::Float64
    total_rate_switch::Float64
    
    # --- PARTICLE COUNTERS (Snapshots) ---
    count_active::Int  # Total State 1 (Mobile + Jammed)
    count_paused::Int  # Total State 2
    count_mobile::Int  # Subset of Active with rate > 0 (The rest are Jammed)
    
    # --- OBSERVABLES (Time-Averaged) ---
    flux_termination::Int # Total Flux
    
    # 1. Particle Type Integrals (For Microscopic Density)
    cum_active_time::Float64   # Integral of N_active
    cum_paused_time::Float64   # Integral of N_paused
    cum_mobile_time::Float64   # Integral of N_mobile
    
    # 2. System State Integrals (For Macroscopic "Regimes")
    
    # Regime A: "Unpaused State" (System has exactly 0 paused particles)
    cum_time_unpaused::Float64 # How long were we in this state?
    cum_mass_unpaused::Float64 # Sum of (Total N * dt) while in this state
    flux_unpaused::Int         # Flux events that happened while in this state
    
    # Regime B: "Paused State" (System has >= 1 paused particle)
    # Time spent here = (Total Time - cum_time_unpaused)
    cum_mass_paused::Float64   # Sum of (Total N * dt) while in this state
    flux_paused::Int           # Flux events that happened while in this state
end

function SimState(model::TranscriptModel)
    return SimState(
        0.0, 0,
        zeros(Int, model.L),
        zeros(UInt8, model.L),
        zeros(Float64, model.L),
        0.0,
        zeros(Float64, model.L),
        0.0, 0.0,
        
        0, 0, 0, # Counts
        
        0,       # Total Flux
        0.0, 0.0, 0.0, # Particle Integrals
        
        0.0, 0.0, 0,   # Unpaused System Stats
        0.0, 0         # Paused System Stats
    )
end