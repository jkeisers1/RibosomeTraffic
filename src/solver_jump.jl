using JumpProcesses
using SciMLBase

"""
    make_jump_problem(model::TranscriptModel)

Creates a SciML JumpProblem compliant with the DifferentialEquations.jl ecosystem.
"""
function make_jump_problem(model::TranscriptModel)
    # 1. Setup
    L = model.L
    l_rib = model.l_ribosome

    # Initial State: Vector of Integers (0=Empty, 1=Active, 2=Paused)
    u0 = zeros(Int, L)
    tspan = (0.0, 100.0) # Default, user overrides this later

    # We collect all physical events as "Jumps"
    jumps = Vector{ConstantRateJump}()

    # ==========================================================================
    # A. INITIATION (Enter at site 1)
    # ==========================================================================
    rate_init = (u, p, t) -> begin
        # Scan the "footprint" zone. If ANY site 1..l_rib is occupied, we are blocked.
        # Uses @view to avoid allocation.
        is_blocked = any(x -> x > 0, @view u[1:p.l_ribosome])
        return is_blocked ? 0.0 : p.α
    end

    affect_init!(integrator) = (integrator.u[1] = 1)
    push!(jumps, ConstantRateJump(rate_init, affect_init!))

    # ==========================================================================
    # B. ELONGATION (Move i -> i+1)
    # ==========================================================================
    for i = 1:(L-1)
        # CLOSURE: We use 'let' to explicitly capture 'i' for this specific jump
        let site = i, next = i + 1, check_idx = i + l_rib

            rate_elong = (u, p, t) -> begin
                # 1. Source: Must have an Active ribosome (1)
                if u[site] != 1
                    return 0.0
                end

                # 2. Check exclusion
                # If blocked by particle l sites ahead no moving.
                if check_idx <= p.L && u[check_idx] != 0
                    return 0.0
                end

                return p.k_elong[site]
            end

            affect_elong!(integrator) =
                (integrator.u[site] = 0; integrator.u[next] = 1; nothing)
            push!(jumps, ConstantRateJump(rate_elong, affect_elong!))
        end
    end

    # ==========================================================================
    # C. TERMINATION (Leave at L)
    # ==========================================================================
    rate_term = (u, p, t) -> (u[L] == 1) ? p.β : 0.0
    affect_term!(integrator) = (integrator.u[L] = 0; nothing)
    push!(jumps, ConstantRateJump(rate_term, affect_term!))

    # ==========================================================================
    # D. PAUSING & UNPAUSING (State Switching)
    # ==========================================================================
    # NOTE: We add these jumps even if k_pause is 0.0, so we can support
    # parameter sweeping (changing k_pause > 0 later) without rebuilding the problem.
    for i = 1:L
        let site = i
            # Active (1) -> Paused (2)
            rate_pause = (u, p, t) -> (u[site] == 1) ? p.k_pause : 0.0
            affect_pause!(integrator) = (integrator.u[site] = 2; nothing)
            push!(jumps, ConstantRateJump(rate_pause, affect_pause!))

            # Paused (2) -> Active (1)
            rate_recov = (u, p, t) -> (u[site] == 2) ? p.k_unpause : 0.0
            affect_recov!(integrator) = (integrator.u[site] = 1; nothing)
            push!(jumps, ConstantRateJump(rate_recov, affect_recov!))
        end
    end

    # ==========================================================================
    # 3. ASSEMBLY
    # ==========================================================================
    # Define the discrete problem
    dprob = DiscreteProblem(u0, tspan, model)

    # Use Gillespie's Direct Method
    # Note: We splat the `jumps` vector into the constructor
    return JumpProblem(dprob, Direct(), jumps...)
end

"""
    run_sciml_simulation(model, t_max)

Helper to run the simulation using the standard SciML pipeline.
"""
function run_sciml_simulation(model::TranscriptModel, t_max::Float64)
    prob = make_jump_problem(model)
    prob = remake(prob, tspan = (0.0, t_max))
    return solve(prob, SSAStepper())
end
