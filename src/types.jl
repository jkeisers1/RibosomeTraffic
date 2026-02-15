
# types for the ribosome model to stop the if statements (multiple Dispatch! :))) )
abstract type RibosomeState end
struct EmptySite <: RibosomeState end
struct Active <: RibosomeState end
struct Paused <: RibosomeState end
struct Jammed <: RibosomeState end

struct TranscriptModel{T<:Real,V<:AbstractVector{T}}
    L::Int
    l_ribosome::Int
    α::T              # Initiation rate
    β::T              # Termination rate
    k_elong::V
    k_pause::T
    k_unpause::T
    delta::T          # Degradation rate
end

# This constructors makes sure that the compiler specifies the types because T could be many as T is Real
function TranscriptModel(L, l_rib, α, β, k_elong, k_pause, k_unpause, delta)
    # Promote all scalar inputs to the same type T (usually Float64)
    T = promote_type(typeof(α), typeof(β), eltype(k_elong), typeof(delta))
    # Convert fields to T
    return TranscriptModel{T,Vector{T}}(
        L,
        l_rib,
        T(α),
        T(β),
        Vector{T}(k_elong),
        T(k_pause),
        T(k_unpause),
        T(delta),
    )
end

mutable struct SimState{T<:Real}
    time::T
    step_count::Int

    lattice::Vector{Int}
    internal_states::Vector{UInt8} # 0=Empty, 1=Active, 2=Paused

    rate_elong::Vector{T}
    rate_initiation::T
    rate_switch::Vector{T}

    total_rate_elong::T
    total_rate_switch::T

    # count active paused and mobile ribosomes
    count_active::Int
    count_paused::Int
    count_mobile::Int

    flux_termination::Int

    # how much active time the systems spends in the state
    cum_active_time::T
    cum_paused_time::T
    cum_mobile_time::T

    # how much time and mass the system transported in the unpaused and paused states
    cum_time_unpaused::T
    cum_mass_unpaused::T
    flux_unpaused::Int

    cum_mass_paused::T
    flux_paused::Int
end


function SimState(model::TranscriptModel{T}) where {T}
    return SimState{T}(
        zero(T),
        0,
        zeros(Int, model.L),
        zeros(UInt8, model.L),
        zeros(T, model.L),
        zero(T),
        zeros(T, model.L),
        zero(T),
        zero(T),
        0,
        0,
        0,
        0,
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        0,
        zero(T),
        0,
    )
end

function reset_statistics!(state::SimState{T}) where {T}
    state.time = zero(T)
    state.step_count = 0

    # Reset all counters
    state.flux_termination = 0
    state.flux_unpaused = 0
    state.flux_paused = 0

    # Reset all integrals
    state.cum_active_time = zero(T)
    state.cum_paused_time = zero(T)
    state.cum_mobile_time = zero(T)
    state.cum_time_unpaused = zero(T)
    state.cum_mass_unpaused = zero(T)
    state.cum_mass_paused = zero(T)

    @info "Statistics reset. Ready for production run."
end
