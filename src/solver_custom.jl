using Random
using Distributions

# --- 1. Rate Updates (Robust Jamming Logic) ---
function update_local_rates!(state::SimState, model::TranscriptModel, i::Int)
    if i < 1 || i > model.L; return; end

    # Check OLD mobile status
    was_mobile = (state.rate_elong[i] > 0.0)

    state.total_rate_elong  -= state.rate_elong[i]
    state.total_rate_switch -= state.rate_switch[i]
    
    if state.lattice[i] == 0
        state.rate_elong[i] = 0.0; state.rate_switch[i] = 0.0
        if was_mobile; state.count_mobile -= 1; end
        return
    end

    status = state.internal_states[i]
    
    # Switching
    current_switch_rate = (status == 1) ? model.k_pause : ((status == 2) ? model.k_unpause : 0.0)
    state.rate_switch[i] = current_switch_rate
    state.total_rate_switch += current_switch_rate

    # Elongation
    current_elong_rate = 0.0
    if status == 1 # Active
        if i == model.L
            current_elong_rate = model.β
        else
            check_index = i + model.l_ribosome
            is_blocked = (check_index <= model.L && state.lattice[check_index] == 1)
            if !is_blocked; current_elong_rate = model.k_elong[i]; end
        end
    end
    state.rate_elong[i] = current_elong_rate
    state.total_rate_elong += current_elong_rate
    
    # Check NEW mobile status
    is_mobile = (current_elong_rate > 0.0)
    
    if was_mobile && !is_mobile
        state.count_mobile -= 1
    elseif !was_mobile && is_mobile
        state.count_mobile += 1
    end
end

function update_initiation_rate!(state::SimState, model::TranscriptModel)
    is_blocked = false
    for i in 1:model.l_ribosome
        if i <= model.L && state.lattice[i] == 1
             is_blocked = true; break
        end
    end
    state.rate_initiation = is_blocked ? 0.0 : model.α
end

# --- 2. Events (Splitting Flux) ---
function move_ribosome!(state::SimState, model::TranscriptModel, i::Int)
    state.lattice[i] = 0
    current_status = state.internal_states[i]
    state.internal_states[i] = 0
    next_site = i + 1
    
    if next_site > model.L
        # Termination
        state.flux_termination += 1
        state.count_active -= 1 
        
        # SPLIT FLUX: Was the system clean or dirty when this happened?
        if state.count_paused == 0
            state.flux_unpaused += 1
        else
            state.flux_paused += 1
        end
        
        update_local_rates!(state, model, i)
        update_local_rates!(state, model, i - model.l_ribosome)
    else
        # Translocation
        state.lattice[next_site] = 1
        state.internal_states[next_site] = current_status
        update_local_rates!(state, model, i)
        update_local_rates!(state, model, next_site)
        update_local_rates!(state, model, i - model.l_ribosome)
        update_local_rates!(state, model, next_site - model.l_ribosome)
        if i <= model.l_ribosome + 1; update_initiation_rate!(state, model); end
    end
end

function switch_ribosome_state!(state::SimState, model::TranscriptModel, i::Int)
    current = state.internal_states[i]
    if current == 1
        state.internal_states[i] = 2
        state.count_active -= 1; state.count_paused += 1
    elseif current == 2
        state.internal_states[i] = 1
        state.count_paused -= 1; state.count_active += 1
    end
    update_local_rates!(state, model, i)
end

# --- 3. Main Loop (Splitting Density) ---
function step!(state::SimState, model::TranscriptModel)
    a0 = state.rate_initiation + state.total_rate_elong + state.total_rate_switch + model.delta
    if a0 <= 1e-16; return; end
    dt = rand(Exponential(1.0 / a0))

    # A. Particle Type Integrals (The "What" are they?)
    state.cum_active_time += Float64(state.count_active) * dt
    state.cum_paused_time += Float64(state.count_paused) * dt
    state.cum_mobile_time += Float64(state.count_mobile) * dt

    # B. System State Integrals (The "Where" are we?)
    # We count TOTAL particles (Active + Paused) in the current regime
    current_total_N = Float64(state.count_active + state.count_paused)
    
    if state.count_paused == 0
        # Unpaused Regime
        state.cum_time_unpaused += dt
        state.cum_mass_unpaused += current_total_N * dt
    else
        # Paused Regime
        state.cum_mass_paused += current_total_N * dt
    end

    state.time += dt
    state.step_count += 1
    
    r = rand() * a0
    
    if r < model.delta
        fill!(state.lattice, 0); fill!(state.internal_states, 0)
        state.count_active = 0; state.count_paused = 0; state.count_mobile = 0
        fill!(state.rate_elong, 0.0); fill!(state.rate_switch, 0.0)
        state.total_rate_elong = 0.0; state.total_rate_switch = 0.0
        state.rate_initiation = model.α
        
    elseif r < model.delta + state.rate_initiation
        state.lattice[1] = 1; state.internal_states[1] = 1
        state.count_active += 1; state.rate_initiation = 0.0
        update_local_rates!(state, model, 1)
        
    elseif r < model.delta + state.rate_initiation + state.total_rate_elong
        r_elong = r - (model.delta + state.rate_initiation)
        cumulative = 0.0; selected_site = 0
        for i in 1:model.L
            rate = state.rate_elong[i]
            if rate > 0
                cumulative += rate
                if cumulative >= r_elong; selected_site = i; break; end
            end
        end
        if selected_site > 0; move_ribosome!(state, model, selected_site); end
        
    else
        r_switch = r - (model.delta + state.rate_initiation + state.total_rate_elong)
        cumulative = 0.0; selected_site = 0
        for i in 1:model.L
            rate = state.rate_switch[i]
            if rate > 0
                cumulative += rate
                if cumulative >= r_switch; selected_site = i; break; end
            end
        end
        if selected_site > 0; switch_ribosome_state!(state, model, selected_site); end
    end
end

function run_custom_simulation(model::TranscriptModel, t_max::Float64)
    state = SimState(model)
    update_initiation_rate!(state, model)
    while state.time < t_max
        step!(state, model)
    end
    return state
end