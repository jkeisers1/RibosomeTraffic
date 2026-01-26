# src/solver_custom.jl
using Random
using Distributions

"""
    update_local_rates!(state, model, i)

Efficiently recalculates the elongation and switching rates for a ribosome at site `i`.
This replaces the expensive `check_jammed` scans.
"""
function update_local_rates!(state::SimState, model::TranscriptModel, i::Int)
    if i < 1 || i > model.L; return; end

    # Reset old rates
    state.total_rate_elong  -= state.rate_elong[i]
    state.total_rate_switch -= state.rate_switch[i]
    
    if state.lattice[i] == 0
        state.rate_elong[i] = 0.0
        state.rate_switch[i] = 0.0
        return
    end

    status = state.internal_states[i]
    
    # 1. Switching Rates (Same as before)
    current_switch_rate = 0.0
    if status == 1; current_switch_rate = model.k_pause;
    elseif status == 2; current_switch_rate = model.k_unpause; end
    state.rate_switch[i] = current_switch_rate
    state.total_rate_switch += current_switch_rate

    # 2. Elongation / Termination Rates
    current_elong_rate = 0.0
    
    if status == 1 # Only mobile particles move
        if i == model.L
            # --- TERMINATION CASE ---
            # At the last site, we don't check exclusion. We just leave.
            current_elong_rate = model.β
        else
            # --- BULK ELONGATION CASE ---
            # Check footprint ahead
            check_index = i + model.l_ribosome
            is_blocked = false
            
            # If the "nose" of the ribosome is still on the lattice, check for collision
            if check_index <= model.L && state.lattice[check_index] == 1
                is_blocked = true
            end
            
            if !is_blocked
                current_elong_rate = model.k_elong[i]
            end
        end
    end
    
    state.rate_elong[i] = current_elong_rate
    state.total_rate_elong += current_elong_rate
end
"""
    move_ribosome!(state, model, i)

Executes a move from `i` to `i+1`. Handles termination if `i == L`.
Crucially, it triggers rate updates for neighbors.
"""
function move_ribosome!(state::SimState, model::TranscriptModel, i::Int)
    # 1. Pick up the particle
    state.lattice[i] = 0
    current_status = state.internal_states[i]
    state.internal_states[i] = 0
    
    next_site = i + 1
    
    if next_site > model.L
        # ====================================================
        # TERMINATION EVENT (Exit)
        # ====================================================
        state.flux_termination += 1
        state.count_active -= 1 
        
        # 1. The last site `i` is now empty. Set its rate to 0.
        update_local_rates!(state, model, i)
        
        # 2. The ribosome BEHIND us (at L - l_ribosome) was blocked.
        #    Now that we are gone, it sees empty space. Wake it up!
        prev_site = i - model.l_ribosome
        update_local_rates!(state, model, prev_site)
        
    else
        # ====================================================
        # TRANSLOCATION EVENT (Move i -> i+1)
        # ====================================================
        state.lattice[next_site] = 1
        state.internal_states[next_site] = current_status
        
        # Update rates at OLD site (empty), NEW site (occupied), and BEHIND (unblocked)
        update_local_rates!(state, model, i)
        update_local_rates!(state, model, next_site)
        update_local_rates!(state, model, i - model.l_ribosome)
        
        # NEW: Check if we just blocked someone at the new position
        update_local_rates!(state, model, next_site - model.l_ribosome)
        
        # NEW: Check if we unblocked the initiation site
        if i <= model.l_ribosome + 1
            update_initiation_rate!(state, model)
        end
    end
end
"""
    update_initiation_rate!(state, model)

Checks if the first `l_ribosome` sites are empty. 
If yes, initiation is possible (rate = α). If no, rate = 0.
"""
function update_initiation_rate!(state::SimState, model::TranscriptModel)
    # Check if the "landing pad" is clear
    is_blocked = false
    
    # We check sites 1 to l_ribosome
    for i in 1:model.l_ribosome
        if i <= model.L && state.lattice[i] == 1 # There is a ribosome head here
             is_blocked = true
             break
        end
        # Note: In a more complex model where lattice[i] stores 
        # occupancy of *any* part of the ribosome, we'd sum them. 
        # Since your lattice stores "heads" (1) and we assume body trails/leads,
        # we strictly need to ensure no other ribosome head is in [1, l].
        # For simplicity, let's assume `lattice[i] == 1` means a HEAD is at i.
    end
    
    # Update the rate
    state.rate_initiation = is_blocked ? 0.0 : model.α
end
"""
    step!(state, model)

The main Gillespie loop step. Finds the next event and executes it.
"""
function step!(state::SimState, model::TranscriptModel)
    # 1. Total Propensity
    a0 = state.total_rate_elong + state.total_rate_switch + state.rate_initiation + model.delta
    
    # 2. Time Step
    if a0 <= 1e-16 # Float safety for 0
        return # Deadlock
    end
    dt = rand(Exponential(1.0 / a0))

    # Add (Number of Particles * Duration) to the accumulator
    state.cum_active_time += Float64(state.count_active) * dt
    state.cum_paused_time += Float64(state.count_paused) * dt

    # all particles that can move contribute to total elongation rate
    state.cum_moving_masses += state.total_rate_elong * dt

    state.time += dt
    state.step_count += 1
    
    # 3. Select Reaction Channel
    r = rand() * a0
    # --- BRANCH 1: INITIATION ---
    if r < model.delta
        # 1. Wipe the lattice (Ribosomes fall off / Abort)
        fill!(state.lattice, 0)
        fill!(state.internal_states, 0)
        
        # 2. Reset Counters (They are all gone!)
        state.count_active = 0
        state.count_paused = 0
        
        # 3. Reset Rates (Lattice is empty, so no elongation/switching)
        fill!(state.rate_elong, 0.0)
        fill!(state.rate_switch, 0.0)
        state.total_rate_elong = 0.0
        state.total_rate_switch = 0.0
        
        # 4. Immediate Rebirth (New mRNA is ready)
        # The lattice is empty, so initiation is unblocked.
        state.rate_initiation = model.α

    elseif r < state.rate_initiation + model.delta
        # Execute Initiation
        # Place a new mobile ribosome at site 1
        state.lattice[1] = 1
        state.internal_states[1] = 1 # Mobile
        state.count_active += 1
        
        # The entrance is now blocked!
        state.rate_initiation = 0.0
        
        # Calculate rates for this new particle
        update_local_rates!(state, model, 1)
        
    # --- BRANCH 2: ELONGATION ---
    elseif r < state.rate_initiation + state.total_rate_elong + model.delta
        # Shift r to be relative to the elongation block
        r_elong = r - (state.rate_initiation + model.delta)
        cumulative = 0.0
        selected_site = 0
        for i in 1:model.L
            rate = state.rate_elong[i]
            if rate > 0
                cumulative += rate
                if cumulative >= r_elong
                    selected_site = i
                    break
                end
            end
        end
        move_ribosome!(state, model, selected_site)
        
    else
        # --- SWITCHING ---
        r_switch = r - (state.total_rate_elong + state.rate_initiation + model.delta)
        cumulative = 0.0
        selected_site = 0
        for i in 1:model.L
            rate = state.rate_switch[i]
            if rate > 0
                cumulative += rate
                if cumulative >= r_switch
                    selected_site = i
                    break
                end
            end
        end
        switch_ribosome_state!(state, model, selected_site)
    end
end

function switch_ribosome_state!(state::SimState, model::TranscriptModel, i::Int)
    current = state.internal_states[i]
    
    # Flip the state
    if current == 1 # Mobile -> Paused
        state.internal_states[i] = 2
        state.count_active -= 1
        state.count_paused += 1
    elseif current == 2 # Paused -> Mobile
        state.internal_states[i] = 1
        state.count_paused -= 1
        state.count_active += 1
    end
    
    # Recalculate rates for this particle
    update_local_rates!(state, model, i)
end

"""
    run_custom_simulation(model::TranscriptModel, t_max::Float64)

Runs the simulation starting from an empty lattice until `t_max`.
Returns the final `SimState`.
"""
function run_custom_simulation(model::TranscriptModel, t_max::Float64)
    # 1. Initialize empty state
    state = SimState(model)
    
    # 2. CRITICAL: Calculate initial rates
    # The lattice is empty, so initiation is possible.
    update_initiation_rate!(state, model)
    
    # 3. The Loop
    while state.time < t_max
        step!(state, model)
    end
    
    return state
end