# src/types.jl
using StaticArrays

"""
    TranscriptModel

Immutable struct holding all static parameters of the simulation.
"""
struct TranscriptModel
    L::Int                  # Length of lattice (codons)
    l_ribosome::Int         # Footprint size (width in codons, e.g., 10)
    
    # Rates
    α::Float64              # Initiation rate (at site 1)
    β::Float64              # Termination rate (at site L)
    k_elong::Vector{Float64}# Elongation rates per site (vector of length L)
    
    # Pausing Kinetics
    k_pause::Float64        # Rate Mobile -> Paused (k₋)
    k_unpause::Float64      # Rate Paused -> Mobile (k₊)
end

"""
    SimState

Mutable struct holding the changing state of the simulation.
Designed for zero-allocation updates.
"""
mutable struct SimState
    time::Float64
    step_count::Int
    
    # Lattice: 0 = empty, 1 = ribosome head position. 
    # The ribosome body occupies [i, i + l - 1].
    lattice::Vector{Int} 
    
    # Internal State: 0 = None, 1 = Mobile, 2 = Paused
    internal_states::Vector{UInt8} 
    
    # Propensities (Rates)
    # rate_elong[i]: Rate of particle at i moving to i+1
    rate_elong::Vector{Float64}
    
    rate_initiation::Float64   # Rate of initiation at site 1

    # rate_switch[i]: Rate of particle at i switching state (Mobile <-> Paused)
    rate_switch::Vector{Float64}
    
    # Total Propensities (Sums maintained incrementally)
    total_rate_elong::Float64
    total_rate_switch::Float64
    
    # Observables / Counters
    flux_termination::Int
    count_active::Int
    count_paused::Int

    #Observables (Time-averaged Integrals)
    cum_active_time::Float64
    cum_paused_time::Float64
    cum_moving_masses::Float64
end

# Helper constructor to initialize state easily
# Helper constructor to initialize state easily
function SimState(model::TranscriptModel)
    return SimState(
        0.0,                    # time
        0,                      # step_count
        zeros(Int, model.L),    # lattice
        zeros(UInt8, model.L),  # internal_states
        zeros(Float64, model.L),# rate_elong
        0.0,                    # rate_initiation (FIXED: Added this missing 0.0)
        zeros(Float64, model.L),# rate_switch
        0.0,                    # total_rate_elong
        0.0,                    # total_rate_switch
        0,                      # flux_termination
        0,                      # count_mobile
        0,                       # count_paused
        0.0,                    # cum_active_time
        0.0,                    # cum_paused_time
        0.0                     # cum_moving_masses
    )
end