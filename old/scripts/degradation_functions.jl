using SpecialFunctions # Required for erf

"""
    calculate_T0(alpha, epsilon, L, kp, tau; ell=10)

Calculates the expected time an mRNA spends in the unpaused state (T0),
accounting for the competition between stochastic pausing and mRNA degradation.

# Arguments
- `alpha`: Initiation rate
- `epsilon`: Elongation rate
- `L`: Lattice length (codons)
- `kp`: Pausing rate (antibiotic concentration dependent)
- `tau`: Mean mRNA lifetime
- `ell`: Ribosome footprint (default 10)
"""
function calculate_T0(alpha, epsilon, L, kp, tau; ell=10)
    # 1. Calculate Standard TASEP Parameters
    v = epsilon - alpha
    # Time to reach steady state (domain wall travel time)
    tL = L / v 
    
    # Steady state current and density (unpaused)
    J0 = (alpha * (epsilon - alpha)) / (epsilon + alpha * (ell - 1))
    rho0 = alpha / (epsilon + alpha * (ell - 1))

    # 2. Define Integration Constants
    # From Eq. 10: Lambda(t) = a*t^2 for t < tL
    a = (J0 * kp) / 2
    # From Eq. 10: Lambda(t) linear slope b for t >= tL
    b = rho0 * L * kp
    
    # 3. Calculate Part 1: Integral from 0 to tL (Shifted Gaussian)
    # We solve integral of exp(-a*t^2 - t/tau)
    # Using the identity involving erf and completing the square
    
    # Pre-factor from completing the square
    pre_factor_exp = exp(1 / (4 * a * tau^2))
    factor_front = (sqrt(pi) / (2 * sqrt(a))) * pre_factor_exp
    
    # Error function arguments
    u_lower = 1 / (2 * tau * sqrt(a))
    u_upper = sqrt(a) * tL + u_lower
    
    part1 = factor_front * (erf(u_upper) - erf(u_lower))

    # 4. Calculate Part 2: Integral from tL to infinity (Exponential Decay)
    # We solve integral of exp(-a*tL^2 - b(t-tL) - t/tau)
    # This simplifies to a simple exponential integral
    
    decay_rate = b + (1 / tau)
    numerator = exp(-a * tL^2 - (tL / tau))
    
    part2 = numerator / decay_rate

    # 5. Total Effective Unpaused Time
    return part1 + part2
end

"""
    calculate_Tp_total(L, kp, km, l_ribosome, init, elong, tau)

Calculates the expected total time an mRNA spends in the paused state (Tp)
accounting for:
1. The probability that the mRNA degrades BEFORE ever pausing.
2. The race between cluster dissolution and mRNA degradation once paused.

Formula: Tp = P(entry) * Tp_duration_effective
"""
function calculate_Tp_total(L, kp, km, l_ribosome, init, elong, tau)
    # 1. Get Effective Unpaused Time (T0)
    # This must calculate the integral of S(t)*exp(-t/tau)
    T0_eff = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)

    # 2. Probability of Entering the Paused State
    # If we spent T0_eff unpaused out of a total life tau, the rest is "lost" to the pause
    # (or rather, the pause happened before death)
    P_entry = 1.0 - (T0_eff / tau)

    # If probability is effectively zero, return 0 to avoid NaNs
    if P_entry <= 0.0 || kp == 0.0
        return 0.0
    end

    # 3. Calculate Properties of the Jam
    # Get the average position where the first pause occurs
    avg_pos = analyitcal_total_expected_pausing_position_extended(init, elong, L, kp, l_ribosome)
    
    # Initial cluster size (particles)
    Ni = avg_pos / l_ribosome
    
    # Average Batch Size (B)
    # Using your existing function 'distance_paused_particles'
    d = distance_paused_particles(km, kp, Ni)
    B = d + 1.0
    
    # Number of batches required to clear (if life was infinite)
    n_req = Ni / B
    
    # Time to clear (if life was infinite)
    T_clear_inf = n_req / km 

    # 4. Effective Duration of the Jam
    # The expected duration limited by exponential decay tau
    Tp_dur = tau * (1.0 - exp(-T_clear_inf / tau))

    # 5. Total Expected Paused Time
    return P_entry * tau#* Tp_dur

end

"""
    calculate_avg_pause_position_robust(alpha, epsilon, L, kp, tau; ell=10)

Calculates the expected First Pausing Position <X_pause>.
Includes a stability check for low kp to prevent numerical overflow.
"""
function calculate_avg_pause_position_robust(alpha, epsilon, L, kp, tau; ell=10)
    # Edge case: No antibiotics or instant death
    if kp <= 1e-20 || tau <= 0.0
        return L / 2.0 
    end

    # --- 1. Physics Parameters ---
    v = epsilon - alpha
    if v <= 1e-9; v = epsilon; end
    tL = L / v
    
    # Calculate Steady State Density rho0
    if alpha/epsilon >= 1/ (1 + sqrt(ell))
        rho0 = 1 / (ell + sqrt(ell))
    else
        rho0 = alpha / (epsilon + alpha * (ell - 1))
    end
    J0 = rho0 * v

    # Integration constants
    # Rate during filling: k(t) = 2at
    a = (J0 * kp) / 2.0
    c = 1.0 / tau
    b = rho0 * L * kp

    # --- 2. STABILITY CHECK (The Fix) ---
    # We check the exponent term: 1 / (4*a*tau^2)
    # If this is > 100, exp() > 2e43, which causes precision loss with erf.
    # We switch to the "Linear Limit" (a -> 0) approximation.
    
    exponent_check = 1.0 / (4.0 * a * tau^2)
    
    if exponent_check > 10.0
        # === LOW KP BRANCH (Stable) ===
        # Here we assume exp(-at^2) approx 1.0 (Gaussian curvature is negligible)
        # We only integrate against the exponential decay exp(-t/tau).
        
        # 1. Weights (Denominator)
        # Weight_Fill = Int_0^tL (2at) * exp(-ct) dt
        # = 2a * [ (-t/c - 1/c^2)*exp(-ct) ]_0^tL
        term_L_w = exp(-c*tL) * (-tL/c - 1/c^2)
        term_0_w = -1/c^2
        W_fill = 2*a * (term_L_w - term_0_w)
        
        # Weight_Steady = b * Int_tL^Inf exp(-ct) dt
        W_steady = b * (exp(-c*tL) / c)
        
        # 2. Weighted Positions (Numerator)
        # Pos_Fill = Int_0^tL (vt/2) * (2at) * exp(-ct) dt
        # = a*v * Int t^2 * exp(-ct) dt
        # Int t^2 exp(-ct) = [ exp(-ct)/(-c) * (t^2 + 2t/c + 2/c^2) ]
        # Note: Integral x^2 e^kx = e^kx/k (x^2 - 2x/k + 2/k^2). Here k = -c.
        
        inv_neg_c = -1.0/c
        term_L_p = exp(-c*tL) * inv_neg_c * (tL^2 + 2*tL/c + 2/c^2)
        term_0_p = inv_neg_c * (0 + 0 + 2/c^2)
        
        Pos_fill = (a * v) * (term_L_p - term_0_p)
        
        # Pos_Steady = (L/2) * Weight_Steady
        Pos_steady = (L / 2.0) * W_steady
        
        return (Pos_fill + Pos_steady) / (W_fill + W_steady)
    end

    # === STANDARD BRANCH (Your original code) ===
    # Only runs when kp is large enough to be stable
    
    shift = c / (2 * a)
    pre_factor_gauss = (sqrt(pi) / (2 * sqrt(a))) * exp(a * shift^2)
    
    u_0 = shift * sqrt(a)
    u_L = (tL + shift) * sqrt(a)
    
    integ_S_fill = pre_factor_gauss * (erf(u_L) - erf(u_0))
    prob_survive_filling = exp(-a * tL^2 - c * tL)
    decay_steady = b + c
    prob_pause_steady = prob_survive_filling * (b / decay_steady)
    
    term_boundary = 1.0 - prob_survive_filling
    prob_pause_filling = term_boundary - (c * integ_S_fill)
    P_entry = prob_pause_filling + prob_pause_steady
    
    if P_entry <= 1e-20; return L/2.0; end

    # Numerator Integrals
    pre_factor_t2 = exp(a * shift^2) / sqrt(a)
    F0_L = (sqrt(pi)/2) * erf(u_L); F0_0 = (sqrt(pi)/2) * erf(u_0)
    F1_L = -0.5 * exp(-u_L^2);      F1_0 = -0.5 * exp(-u_0^2)
    val_F0 = F0_L - F0_0
    val_F1 = F1_L - F1_0
    val_F2 = 0.5 * (val_F0 - (u_L*exp(-u_L^2) - u_0*exp(-u_0^2)))
    
    term_y2 = (1/a) * val_F2
    term_y1 = (-2*shift/sqrt(a)) * val_F1
    term_y0 = (shift^2) * val_F0
    
    integral_t2 = pre_factor_t2 * (term_y2 + term_y1 + term_y0)
    
    pos_filling_contribution = (a * v) * integral_t2
    pos_steady_contribution = (L / 2.0) * prob_pause_steady
    
    return (pos_filling_contribution + pos_steady_contribution) / P_entry
end

"""
    calculate_state_currents(L, kp, km, l_ribosome, init, elong, tau)

Calculates P0, J0, Pp, Jp for the Fold Change model.
- Jp is PURE RUN-OFF (Leakage = 0).
- Includes "Gap Logic" to handle early pauses on long genes.
"""
function calculate_state_currents(L, kp, km, l_ribosome, init, elong, tau)
    # --- 1. Calculate Times & Weights ---
    # T0: Expected time unpaused (weighted by death)
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    
    # Tp_total: Expected time paused (weighted by death + P_entry)
    Tp_total = calculate_Tp_total(L, kp, km, l_ribosome, init, elong, tau)
    
    total_life_calc = T0 + Tp_total
    
    if total_life_calc <= 1e-20
        return 0.0, 0.0, 0.0, 0.0
    end

    # Probabilities (Time Fractions)
    P0 = T0 / total_life_calc
    Pp = Tp_total / total_life_calc
    
    # P_entry: Probability of ever entering the paused state
    P_entry = 1.0 - (T0 / tau)
    if P_entry < 0; P_entry = 0.0; end


    # --- 2. Physics Parameters ---
    v = elong - init; if v <= 1e-9; v = elong; end
    tL = L / v
    
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        J_ss = elong / (1 + sqrt(l_ribosome))^2
        rho_ss = 1 / (l_ribosome + sqrt(l_ribosome))
    else
        J_ss = (init * (elong - init)) / (elong + init * (l_ribosome - 1))
        rho_ss = init / (elong + init * (l_ribosome - 1))
    end


    # --- 3. Calculate J0 (Unpaused Current) ---
    
    # 'a' parameter for the Rayleigh distribution of first pause times
    a = (rho_ss * v * kp) / 2.0
    
    # Probability of surviving the startup ramp without pausing or dying
    P_reach_steady = exp(-a * tL^2 - tL / tau)
    
    # Expected duration of productive steady state
    rate_steady = (rho_ss * L * kp) + (1.0 / tau)
    T_prod = P_reach_steady / rate_steady
    
    particles_unpaused = J_ss * T_prod
    
    if T0 > 1e-20
        J0 = particles_unpaused / T0
    else
        J0 = 0.0
    end


    # --- 4. Calculate Jp (Paused Current) - RUN-OFF ONLY ---
    
    # A. Geometry: Where is the jam?
    avg_pos = calculate_avg_pause_position_robust(init, elong, L, kp, tau; ell=l_ribosome)
    
    # B. The Two Regimes (Gap vs. Full)
    # Heuristic: The stream front is approx 2x the average pause position during filling
    wall_pos = 2.0 * avg_pos
    
    t_gap_start = 0.0
    t_gap_end = 0.0
    
    if wall_pos < L
        # === FILLING PHASE (Gap Exists) ===
        stream_len = avg_pos 
        gap_len = L - wall_pos
        if stream_len < 0; stream_len = 0; end
        
        # Ribosomes must cross gap_len first (delay)
        t_gap_start = gap_len / v
        # Then the stream takes stream_len/v to exit
        t_gap_end   = (gap_len + stream_len) / v
    else
        # === STEADY PHASE (No Gap) ===
        stream_len = L - avg_pos
        if stream_len < 0; stream_len = 0; end
        
        # Ribosomes exit immediately
        t_gap_start = 0.0
        t_gap_end   = stream_len / v
    end
    
    # C. Calculate Run-Off Particles
    # N_runoff = J_ss * Int_{t_start}^{t_end} exp(-t/tau) dt
    N_runoff = J_ss * tau * (exp(-t_gap_start/tau) - exp(-t_gap_end/tau))
    
    # D. Leakage = 0 (Explicitly Ignored)
    N_leak = 0.0
    
    # E. Convert to Conditional Rate Jp
    # We need the effective duration of *one* pause event to normalize the rate.
    Ni = avg_pos / l_ribosome
    B = (km/kp) + 1.0; n_req = Ni / B; T_clear_inf = n_req / km 
    Tp_dur_single = tau * (1.0 - exp(-T_clear_inf / tau))
    
    if Tp_dur_single > 1e-20
        # Jp is the rate of particles exiting *while the system is in state P*.
        # Jp = Total_Exiting_Particles / Duration_of_State
        Jp = (N_runoff + N_leak) / Tp_dur_single
    else
        Jp = 0.0 
    end

    return P0, J0, Pp, Jp
end

function calculate_p_entry(L, kp, l_ribosome, init, elong, tau)
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    P_entry = 1.0 - (T0 / tau)
    if P_entry < 0; P_entry = 0.0; end
    return P_entry
end

function calculate_P0(L, kp, l_ribosome, init, elong, tau)
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    return T0 / tau
end


#########################################
"""
    calc_rho_paused_part_exact(km, kp, l_ribosome, L, init, elong, tau)

Calculates the Density contribution from the PAUSED state.
Includes:
1. The Jam (Weighted decay)
2. The Run-Off (Downstream particles exiting)
"""
function calc_rho_paused_part_exact(km, kp, l_ribosome, L, init, elong, tau)
    # Safety
    if kp == 0.0 || tau == 0.0; return 0.0; end

    # --- 1. Probability of Entry ---
    # We need to know how likely we are to be in this state
    T0_total = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    if T0_total == 0; return 0.0; end
    P_entry = 1.0 - (T0_total / tau)
    if P_entry <= 1e-9; return 0.0; end

    # --- 2. Determine Geometry (Wall & Stream) ---
    # Use the correct expected pausing position (with decay censorship)
    avg_pos = calculate_Xf_integrated(init, elong, L, kp, l_ribosome, tau)
    
    # 3. Calculate "Particles in the Jam" (Your existing logic)
    Ni = avg_pos / l_ribosome
    
    # We use your discrete sum logic to get the time-averaged particles in the jam
    # <N_c>_time_avg = Sum(Particles * Duration) / tau
    # We can reuse the structure of your previous function but adapted for the weighted contribution
    
    # Batch properties
    d = (km / kp) # Assuming Ni is large enough that term goes to 0, or use full formula
    B = d + 1.0
    nB = ceil(Int, Ni / B)
    
    # Survival p
    p = km / (km + 1/tau)
    
    # Weighted Sum of Particles * Time
    # Each batch q lasts for time (1/km)
    # Weight = p^q
    # Particles = (Ni - B*q)
    # We need sum: (1/km) * sum (Ni - Bq) * p^q
    
    if abs(1.0 - p) < 1e-6
        S1 = Float64(nB) # Sum p^q
        S2 = nB*(nB-1.0)/2.0 # Sum q*p^q
    else
        S1 = (1.0 - p^nB) / (1.0 - p)
        term_a = (p * (1.0 - p^(nB - 1))) / ((1.0 - p)^2)
        term_b = ((nB - 1) * p^nB) / (1.0 - p)
        S2 = term_a - term_b
    end
    
    # Integral of Jam Particles over time
    # (Ni * S1 - B * S2) is the sum of particles weighted by survival steps
    # Multiply by (1/km) to get "Particle-Seconds"
    integral_Jam = (Ni * S1 - B * S2) * (1.0 / km)

    # --- 4. Calculate "Run-Off Particles" ---
    # These are the particles downstream of the jam.
    
    v = elong - init
    if init/elong >= 1/(1+sqrt(l_ribosome))
        rho_ss = 1 / (l_ribosome + sqrt(l_ribosome))
    else
        rho_ss = init / (elong + init * (l_ribosome - 1))
    end
    
    # Geometry Check (Filling vs Steady)
    wall_pos = 2.0 * avg_pos
    
    if wall_pos < L
        # Filling Phase: Stream length is avg_pos. Gap is L - wall_pos.
        stream_len = avg_pos
        gap_len = L - wall_pos
        t_start = gap_len / v
        t_end = (gap_len + stream_len) / v
    else
        # Steady Phase: Stream length is L - avg_pos. No gap.
        stream_len = L - avg_pos
        gap_len = 0.0
        t_start = 0.0
        t_end = stream_len / v
    end
    
    # We need Integral of N_runoff(t) * exp(-t/tau) dt
    # N_runoff(t) decreases linearly as particles exit.
    # N(t) = rho_ss * (stream_len - v * (t - t_start))  for t in [t_start, t_end]
    
    # This integral is solvable: Int (A - Bt) * exp(-t/tau)
    # = tau * exp(-t/tau) * (A - B(t + tau))
    
    if stream_len > 0
        A = rho_ss * (stream_len + v * t_start) # Constant part to make math easier
        B = rho_ss * v
        
        # Function to evaluate antiderivative: F(t) = -tau * exp(-t/tau) * (N(t) + B*tau)
        # But N(t) is the count. Let's use the explicit form:
        # Int = [-tau * exp(-t/tau) * (Current_Count(t) + J_ss * tau)] evaluated at limits
        
        J_ss = rho_ss * v
        
        # Value at t_start (Count = rho_ss * stream_len)
        val_start = -tau * exp(-t_start/tau) * ((rho_ss * stream_len) + (J_ss * tau))
        
        # Value at t_end (Count = 0)
        val_end   = -tau * exp(-t_end/tau) * (0.0 + (J_ss * tau))
        
        integral_Runoff = val_end - val_start
    else
        integral_Runoff = 0.0
    end
    
    # --- 5. Final Density Contribution ---
    # Total Particle-Seconds / (L * tau) * P_entry
    # (Dividing by L gives density, Dividing by tau averages over lifetime)
    
    total_integral = integral_Jam + integral_Runoff
    
    # Note: P_entry is the probability we ENTER the paused state.
    # But the integrals above are "Expected Particle-Seconds GIVEN we entered".
    # So we multiply directly.
    
    return P_entry * (total_integral / (L * tau))
end
function calc_J_unpaused_contribution(init, elong, L, kp, l_ribosome, tau)
    if tau == 0.0; return 0.0; end

    # 1. Theoretical Steady State Parameters
    # (Standard TASEP formulas)
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        J_ss = elong / (1 + sqrt(l_ribosome))^2
        rho_ss = 1 / (l_ribosome + sqrt(l_ribosome))
    else
        J_ss = (init * (elong - init)) / (elong + init * (l_ribosome - 1))
        rho_ss = init / (elong + init * (l_ribosome - 1))
    end

    # 2. Define the "Start Line" (t_L)
    v = elong - init
    if v <= 0; v = elong; end # Safety for high initiation
    tL = L / v

    # 3. Probability of Reaching t_L (Surviving the Ramp)
    # The exponent is Integral[0->tL] of (Rate_Pause(t) + Rate_Death)
    # Rate_Pause(t) approx (rho * v * kp * t) -> Integral is a * t^2
    a = (rho_ss * v * kp) / 2.0
    
    P_reach_start = exp(-a * tL^2 - tL / tau)

    # 4. Expected Run Time AFTER t_L
    # Once full, hazard is constant: Pausing (Full Lattice) + Degradation
    rate_steady_decay = (rho_ss * L * kp) + (1.0 / tau)
    
    T_prod = P_reach_start / rate_steady_decay

    # 5. Final Contribution
    # Current * (Productive_Time / Avg_Lifetime)
    return J_ss * (T_prod / tau)
end

"""
    calc_J_paused_part_exact(km, kp, l_ribosome, L, init, elong, tau)

Calculates the Paused Current contribution using the 'Domain Wall' logic.
Handles the 'Gap Delay' that kills current on long genes.
"""
function calc_J_paused_part_exact(km, kp, l_ribosome, L, init, elong, tau)
    # 0. Safety Checks
    if kp == 0.0 || tau == 0.0; return 0.0; end

    # --- 1. Probability of Entering Paused State ---
    T0_total = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    P0 = T0_total / tau
    P_entry = 1.0 - P0
    
    if P_entry <= 1e-9; return 0.0; end

    # --- 2. Determine Geometry (Wall & Gap) ---
    # Where does the pause happen?
    avg_pos = calculate_Xf_integrated(init, elong, L, kp, l_ribosome, tau)
    
    # Speed of the ribosome stream
    v = elong - init 
    
    # Calculate Density using your formula
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        rho = 1 / (l_ribosome + sqrt(l_ribosome))
    else
        rho = init / (elong + init * (l_ribosome - 1))
    end
    
    # Current of the runoff stream
    J_runoff = rho * v

    # --- 3. The Two Regimes ---
    
    # Wall position approximation (Based on your "2 * avg_pos" logic)
    wall_pos = 2.0 * avg_pos
    
    if wall_pos < L
        # === SCENARIO A: FILLING PHASE (Gap Exists) ===
        # The stream is only between avg_pos and wall_pos.
        stream_len = avg_pos 
        
        # The empty gap the first ribosome must cross
        gap_len = L - wall_pos
        
        # Time delays
        t_gap_start = gap_len / v          # Time until first particle exits
        t_gap_end   = (gap_len + stream_len) / v  # Time until last particle exits

    else
        # === SCENARIO B: STEADY STATE (Full Stream) ===
        # The stream goes from avg_pos to L.
        stream_len = L - avg_pos
        
        # No gap
        t_gap_start = 0.0
        t_gap_end   = stream_len / v
    end

    # --- 4. Calculate Expected Particles (Run-Off) ---
    # We integrate the surviving current: Int_{t_start}^{t_end} J * exp(-t/tau) dt
    # This accounts for the mRNA dying while the gap is being crossed.
    
    # Integral of exp(-t/tau) is -tau * exp(-t/tau)
    # Definite integral from t1 to t2: tau * (exp(-t1/tau) - exp(-t2/tau))
    
    particles_runoff = J_runoff * tau * (exp(-t_gap_start/tau) - exp(-t_gap_end/tau))


    # --- 5. Add Trickle (Jam Dissolution) ---
    # (Optional: keeps the small contribution from the jam itself clearing)
    
    # Cluster parameters
    Ni = avg_pos / l_ribosome
    d = distance_paused_particles(km, kp, Ni)
    B = d + 1.0
    nB = ceil(Int, Ni / B)
    
    # Probability of surviving the gap + stream travel
    # (The jam is at avg_pos, so it must travel L - avg_pos to exit)
    dist_jam_exit = L - avg_pos
    t_jam_travel = dist_jam_exit / v
    prob_survive_travel = exp(-t_jam_travel / tau)
    
    # Expected batches cleared
    p = km / (km + 1/tau)
    if abs(1.0 - p) < 1e-6
        expected_batches = Float64(nB)
    else
        expected_batches = p * (1.0 - p^nB) / (1.0 - p)
    end
    
    particles_trickle = (B * expected_batches) * prob_survive_travel


    # --- 6. Total Contribution ---
    total_particles = particles_runoff #+ particles_trickle
    
    # Convert to rate contribution (weighted by entry probability)
    return P_entry * (total_particles / tau)
end


function calculate_density_final(init, l_ribosome, elong, L, kp, km, tau)
    # --- 1. Calculate Weights (Fractions of Lifetime) ---
    T0_eff = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    
    # P0: Fraction of life unpaused
    P0 = T0_eff / tau
    
    # P_entry: Probability of reaching paused state
    P_entry = 1.0 - P0
    if P_entry <= 0; P_entry = 0.0; end

    # Tp_fraction: Fraction of life spent in paused state
    # = P(Entry) * (Duration / tau)
    avg_pos = analyitcal_total_expected_pausing_position_extended(init, elong, L, kp, l_ribosome)
    
    # Duration Calculation
    if kp != 0 
        Ni = avg_pos / l_ribosome
        d = km / kp #distance_paused_particles(km, kp, Ni)
        B = d + 1.0
        nB = ceil(Int, Ni / B)
        T_clear_inf = nB / km
        Tp_dur = tau * (1 - exp(-T_clear_inf / tau))
        
        Pp = P_entry * (Tp_dur / tau)
    else
        Pp = 0.0
    end

    # --- 2. Calculate State Densities ---
    
    # rho_0 (Standard TASEP)
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        rho0 = 1 / (l_ribosome + sqrt(l_ribosome)) 
    else
        rho0 = init / (elong + init * (l_ribosome - 1))
    end
    
    # rho_p (Discrete Weighted Average)
    if kp == 0
        rhoP = 0.0
    else
        rhoP = calculate_discrete_rho_p_closed(km, kp, l_ribosome, avg_pos, tau)
    end

    # --- 3. Weighted Sum ---
    # Matches Paper: rho = P0*rho0 + Pp*rhoP
    return P0 * rho0 + Pp * rhoP
end

"""
    calc_J_unpaused_part(init, l_ribosome, elong, L, kp, tau)

Calculates the current contribution from the Unpaused state.
Formula: J_ss * (Expected Productive Time / tau)
"""
function calc_J_unpaused_part(init, l_ribosome, elong, L, kp, tau)
    # 1. Theoretical J_ss (Infinite Time, Full Lattice)
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        J_ss = elong / (1 + sqrt(l_ribosome))^2
    else
        J_ss = (init * (elong - init)) / (elong + init * (l_ribosome - 1))
    end

    # 2. Calculate Start-Up Time (tL)
    v = elong - init
    tL = L / v

    # 3. Calculate "Productive Time"
    # This is the expected time the mRNA spends alive and unpaused AFTER tL.
    # Integral_{tL}^{inf} S(t) * exp(-t/tau) dt
    
    if kp == 0.0
        # Simple case: S(t) = 1. Integral is just exp decay starting at tL.
        # Int_{tL}^{inf} exp(-t/tau) dt = tau * exp(-tL/tau)
        T_productive = tau * exp(-tL / tau)
    else
        # Complex case: S(t) decays due to pausing.
        # We use the approximation S(t) ≈ exp(-a*t^2) or exp(-b*t) depending on regime.
        # For long genes, the domain wall dominates.
        # Let's use the robust exponential approximation for the probability of reaching steady state.
        
        # 'a' parameter for Rayleigh distribution of pause times
        rho_ss = J_ss / v
        a = (rho_ss * v * kp) / 2.0
        
        # Probability of even reaching the productive phase (surviving tL without pausing or dying)
        P_survive_start = exp(-a * tL^2 - tL / tau)
        
        # Once in steady state, the combined decay rate is (pausing + degradation)
        # rate = (rho_ss * L * kp) + (1/tau)
        b = rho_ss * L * kp
        decay_rate = b + 1/tau
        
        # Expected time in steady state = P(Reach) * (1/rate)
        T_productive = P_survive_start / decay_rate
    end
    
    # 4. Final Contribution
    # Current = J_ss * (Productive_Time / Total_Lifetime)
    if tau == 0.0; return 0.0; end
    
    return J_ss * (T_productive / tau)
end

"""
    calc_J_paused_part(km, kp, l_ribosome, L, init, elong, tau)

Calculates the current contribution from the Paused state.
Formula: P(Entry) * (Expected Particles Cleared / tau)
"""
function calc_J_paused_part(km, kp, l_ribosome, L, init, elong, tau)
    # Safety
    if kp == 0.0 || tau == 0.0; return 0.0; end

    # 1. Calculate P(Entry)
    # Probability that a pause happens before death.
    # P_entry = 1 - (T0_unpaused / tau)
    # Note: We need the TOTAL unpaused time here (including non-productive start-up), 
    # so we use the standard calculate_T0 function.
    T0_total = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    P0 = T0_total / tau
    P_entry = 1.0 - P0
    
    if P_entry <= 1e-9; return 0.0; end

    # 2. Cluster Properties
    avg_pos = analyitcal_total_expected_pausing_position_extended(init, elong, L, kp, l_ribosome)
    Ni = avg_pos / l_ribosome
    d = distance_paused_particles(km, kp, Ni)
    B = d + 1.0
    nB = ceil(Int, Ni / B)
    
    # 3. Calculate Expected Particles Cleared
    # We sum the batches, weighted by the probability of surviving long enough to clear them.
    # p = probability of surviving one batch-time (1/km)
    p = km / (km + 1/tau)
    
    expected_batches = 0.0
    # Sum of Geometric Series: Sum_{q=1}^{nB} p^q
    # = p * (1 - p^nB) / (1 - p)
    if abs(1.0 - p) < 1e-6
        expected_batches = Float64(nB)
    else
        expected_batches = p * (1.0 - p^nB) / (1.0 - p)
    end
    
    expected_particles_cleared = B * expected_batches
    
    # 4. Final Contribution
    # Contribution = P_entry * (Particles / tau)
    # (Dividing by tau effectively spreads these particles over the mRNA lifetime)
    return P_entry * (expected_particles_cleared / tau)
end

function calculate_current_final(init, l_ribosome, elong, L, kp, km, tau)
    # --- 1. Calculate Weights ---
    T0_eff = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    P0 = T0_eff / tau
    
    P_entry = 1.0 - P0
    if P_entry <= 0; P_entry = 0.0; end
    
    # Get Pause Duration for Weight
    avg_pos = analyitcal_total_expected_pausing_position_extended(init, elong, L, kp, l_ribosome)
    Ni = avg_pos / l_ribosome
    d = distance_paused_particles(km, kp, Ni)
    B = d + 1.0
    nB = ceil(Int, Ni / B)
    T_clear_inf = nB / km
    Tp_dur = tau * (1 - exp(-T_clear_inf / tau))
    
    Pp = P_entry * (Tp_dur / tau)

    # --- 2. Calculate State Currents ---
    
    # A. Calculate J0 (Standard TASEP)
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        J0_raw = elong / (1 + sqrt(l_ribosome))^2
    else
        J0_raw = (init * (elong - init)) / (elong + init * (l_ribosome - 1))
    end

    # --- THE FIX FOR LONG LATTICES ---
    # Calculate start-up time tL (Time to first protein)
    v = elong - init
    tL = L / v
    
    # Apply exponential penalty: 
    # Current is only observed if mRNA survives past tL.
    # For L=1350, this reduces J0 by ~40-50%.
    correction = exp(-tL / tau)
    J0 = J0_raw * correction


    # B. Calculate Jp (Effective Paused Current)
    if kp == 0
        Jp = 0.0
    else
        Jp = calculate_discrete_J_p(km, kp, l_ribosome, avg_pos, tau)
    end

    # --- 3. Weighted Sum ---
    return P0 * J0 + Pp * Jp
end

"""
    calculate_density_jam_only_exact(km, kp, l_ribosome, L, init, elong, tau)

Calculates the weighted density contribution of the Paused State based on the 'Jam-Only Model'.
Implements the exact analytical integration of the Jam Growth and Discrete Dissolution phases.

Theory Reference:
- Phase I: The Growing Jam (Linear growth from steady state N_start to N_max).
- Phase II: Discrete Dissolution (Staircase shrinking in batches B).
- Normalization: Weighted by P_entry and total spacetime volume (L * tau).

Returns:
    rho_p_contrib: The value to be added to the unpaused density contribution.
"""
function calculate_density_jam_only_exact2(km, kp, l_ribosome, L, init, elong, tau)
    # Safety checks
    if kp <= 0.0 || tau <= 0.0; return 0.0; end

    # --- 1. Weights & Entry Probability ---
    # We first calculate the time spent unpaused (T0) to determine P_entry.
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    
    # P_entry is the probability that the system effectively enters the paused state 
    # (i.e., a pause happens before the mRNA degrades).
    P_entry = 1.0 - (T0 / tau)
    
    # If probability is effectively zero, return 0.0
    if P_entry <= 1e-20; return 0.0; end


    # --- 2. System Parameters (Steady State) ---
    # Calculate density (rho_ss) and current (J_ss) of the unpaused stream.
    # These determine the initial jam size (N_start) and the filling rate (J_ss).
    
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        # Max Current Phase
        rho_ss = 1 / (l_ribosome + sqrt(l_ribosome))
        J_ss = elong / (1 + sqrt(l_ribosome))^2
    else
        # Low Density Phase
        rho_ss = init / (elong + init * (l_ribosome - 1))
        J_ss = (init * (elong - init)) / (elong + init * (l_ribosome - 1))
    end


    # --- 3. Jam Geometry ---
    # Determine where the jam forms (average pause position X_f).
    avg_pos = calculate_avg_pause_position_robust(init, elong, L, kp, tau; ell=l_ribosome)
    
    # N_start: Number of particles already in the segment when pause starts.
    # Eq: N_start = rho_ss * X_f
    N_start = rho_ss * avg_pos
    
    # N_max: Maximum capacity of the jammed segment.
    # Eq: N_max = X_f / ell
    N_max = avg_pos / l_ribosome


    # --- 4. PHASE I: The Growing Jam ---
    # The jam grows linearly from N_start to N_max at rate J_ss.
    # Eq (1): t_grow = (N_max - N_start) / J_ss
    
    if J_ss > 1e-9
        t_grow = (N_max - N_start) / J_ss
    else
        t_grow = 0.0
    end
    if t_grow < 0; t_grow = 0.0; end 
    
    # E_g: Probability that mRNA survives the growth phase
    E_g = exp(-t_grow / tau)
    
    # Integral of Load during Growth: Int_0^t_grow (N_start + J*t) * exp(-t/tau)
    # Eq (3): L_growth = tau(1 - Eg)N_start + J * [tau^2(1 - Eg) - tau*t_grow*Eg]
    
    term_linear = tau * (1.0 - E_g) * N_start
    term_growth = J_ss * ( (tau^2 * (1.0 - E_g)) - (tau * t_grow * E_g) )
    
    L_growth = term_linear + term_growth


    # --- 5. PHASE II: Discrete Dissolution ---
    # The jam shrinks from N_max down to 0 in steps.
    
    # Batch Size B: How many particles leave per clearing event?
    # Standard TASEP: B = 1. Extended: B = 1 + km/kp.
    d = (km / kp)
    B = d + 1.0 
    
    # Q_max: Max steps to clear
    Q_max = ceil(Int, N_max / B)
    
    # Eq (4): Effective Step Duration (Delta t_eff)
    dt_eff = tau * (1.0 - exp(-1.0 / (km * tau)))
    
    # Eq (below 5): Probability of surviving a single step (p)
    p_step = km / (km + 1.0/tau)
    
    # Calculate Geometric Sums (Eq 6 & 7)
    # Sum1 = Sum_{q=0}^{Q-1} p^q
    # Sum2 = Sum_{q=0}^{Q-1} q * p^q
    
    if abs(1.0 - p_step) < 1e-6
        # Limit p -> 1 (Infinite Life)
        S1 = Float64(Q_max)
        S2 = Q_max * (Q_max - 1.0) / 2.0
    else
        # Standard Finite Life
        pQ = p_step^Q_max
        inv_1_p = 1.0 / (1.0 - p_step)
        
        S1 = (1.0 - pQ) * inv_1_p
        
        term_a = p_step * (1.0 - p_step^(Q_max - 1)) * (inv_1_p^2)
        term_b = (Q_max - 1) * pQ * inv_1_p
        S2 = term_a - term_b
    end
    
    # Eq (5): Total Load during Dissolution (L_diss)
    # L_diss = dt_eff * Sum (N_max - qB) * p^q
    L_diss = dt_eff * (N_max * S1 - B * S2)


    # --- 6. Total Density Contribution ---
    # Eq (8): rho_contrib = (P_entry / L*tau) * (L_growth + E_g * L_diss)
    
    total_load = L_growth + (E_g * L_diss)
    
    rho_p_contrib = (P_entry / (L * tau)) * total_load
    
    return rho_p_contrib
end


function calc_density_paused_total(km, kp, l_ribosome, L, init, elong, tau)
    if kp <= 0.0 || tau <= 0.0; return 0.0; end

    # --- 1. Calculate Weights ---
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    P_entry = 1.0 - (T0 / tau)
    if P_entry <= 1e-20; return 0.0; end

    # --- 2. Part A: The Jam ---
    rho_jam_contrib = calculate_density_jam_only_exact2(km, kp, l_ribosome, L, init, elong, tau)

    # --- 3. Part B: The Run-Off ---
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        rho_ss = 1 / (l_ribosome + sqrt(l_ribosome))
    else
        rho_ss = init / (elong + init * (l_ribosome - 1))
    end
    
    v = elong - init; if v <= 1e-9; v = elong; end

    # Geometry
    avg_pos = calculate_avg_pause_position_robust(init, elong, L, kp, tau; ell=l_ribosome)
    wall_pos = 2.0 * avg_pos
    
    # Define Run-Off Stream
    if wall_pos < L
        # Filling Phase
        stream_len = avg_pos
        gap_len = L - wall_pos
    else
        # Steady Phase
        stream_len = L - avg_pos
        gap_len = 0.0
    end
    
    if stream_len < 0; stream_len = 0.0; end

    t_start_exit = gap_len / v
    t_end_exit   = (gap_len + stream_len) / v
    
    # INTEGRAL 1: Gap Delay
    load_gap = 0.0
    if t_start_exit > 1e-9
        mass_gap = rho_ss * stream_len
        integral_gap = tau * (1.0 - exp(-t_start_exit / tau))
        load_gap = mass_gap * integral_gap
    end
    
    # INTEGRAL 2: Exit Phase
    load_exit = 0.0
    if t_end_exit > t_start_exit
        duration = t_end_exit - t_start_exit
        prob_survive_gap = exp(-t_start_exit / tau)
        
        # Constant Part Integral
        int_const = tau * (1.0 - exp(-duration / tau))
        term_1 = (rho_ss * stream_len) * int_const
        
        # Linear Decay Integral
        c = 1.0/tau
        val_lower = -1.0
        val_upper = exp(-duration/tau) * (-duration/tau - 1.0)
        
        # This term calculates the integral of (u * exp(-u/tau))
        # It is strictly positive.
        term_2 = (rho_ss * v) * (tau^2 * (val_upper - val_lower))
        
        # === THE FIX IS HERE ===
        # Mass = Constant - Linear_Loss
        # Therefore: Integral = Term1 - Term2
        load_exit = prob_survive_gap * (term_1 - term_2)
    end
    
    total_runoff_load = load_gap + load_exit
    
    rho_runoff_contrib = P_entry * (total_runoff_load / (L * tau))

    return rho_jam_contrib + rho_runoff_contrib
end

"""
    calc_rho_unpaused_contribution(init, l_ribosome, elong, L, kp, tau)

Calculates the time-averaged density contribution from the UNPAUSED state (rho_0).
It integrates the density profile N(t)/L over the expected unpaused lifetime.

Phases:
1. Filling Phase (0 to tL): Density grows linearly as the lattice fills.
2. Steady Phase (tL to infinity): Density is constant at rho_ss.

Returns:
    rho_0_contrib: The weighted density value (already multiplied by P0 logic).
"""
function calc_rho_unpaused_contribution(init, l_ribosome, elong, L, kp, tau)
    # Safety Check
    if tau <= 0.0; return 0.0; end

    # --- 1. Steady State Parameters ---
    # Standard TASEP Mean Field Density (rho_ss)
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        rho_ss = 1 / (l_ribosome + sqrt(l_ribosome))
    else
        rho_ss = init / (elong + init * (l_ribosome - 1))
    end

    # --- 2. Timescales ---
    # Velocity of the domain wall (v)
    v = elong - init
    if v <= 1e-9; v = elong; end
    
    # Time to fill the lattice (tL)
    tL = L / v

    # --- 3. Integration Constants ---
    # We are integrating N(t) * S(t) where S(t) is survival probability.
    # Survival S(t) = exp(-a*t^2 - t/tau) during filling.
    
    # 'a' parameter (Rate of pausing during filling)
    # Hazard rate h(t) = rho_ss * v * kp * t
    # Integral H(t) = a * t^2
    a = (rho_ss * v * kp) / 2.0
    
    # 'c' parameter (Degradation rate)
    c = 1.0 / tau


    # --- 4. PHASE I: Filling Phase (0 to tL) ---
    # Density profile: rho(t) = rho_ss * (v*t / L)
    # Total Particles: N(t) = rho_ss * v * t
    # Integral 1: Int_0^tL (rho_ss * v * t) * exp(-a*t^2 - c*t) dt
    
    # This integral is of the form: Int t * exp(-at^2 - ct)
    # We solve this using the derivative of the Gaussian integral or numerical approx.
    
    if a > 1e-9
        # Constants for completing the square
        shift = c / (2 * a)
        sqrt_a = sqrt(a)
        
        # We use the identity: Int t*exp(-at^2 - ct) = (-1/2a) * [exp(-at^2 - ct) + c * Int exp(-at^2 - ct)]
        
        # 1. Boundary term: [exp(-at^2 - ct)]_0^tL
        term_boundary = exp(-a*tL^2 - c*tL) - 1.0
        
        # 2. Gaussian Integral term: Int_0^tL exp(-at^2 - ct)
        # We reuse the logic from the T0 calculation
        u_0 = shift * sqrt_a
        u_L = (tL + shift) * sqrt_a
        pre_factor_gauss = (sqrt(pi) / (2 * sqrt_a)) * exp(a * shift^2)
        integral_gauss = pre_factor_gauss * (erf(u_L) - erf(u_0))
        
        # Combine
        integral_filling = (-1.0 / (2*a)) * (term_boundary + c * integral_gauss)
        
        load_filling = (rho_ss * v) * integral_filling

    else
        # Small 'a' limit (Linear degradation only, ignoring pausing during filling)
        # Int t * exp(-ct) = [exp(-ct)/c^2 * (-ct - 1)]
        val_L = exp(-c * tL) * (-c * tL - 1.0)
        val_0 = -1.0
        integral_filling = (val_L - val_0) / c^2
        
        load_filling = (rho_ss * v) * integral_filling
    end


    # --- 5. PHASE II: Steady State (tL to Infinity) ---
    # Density is constant: rho_ss
    # Total Particles: N_ss = rho_ss * L
    # We just need the expected time spent in steady state.
    
    # Probability of surviving filling phase unpaused
    prob_survive_filling = exp(-a * tL^2 - c * tL)
    
    # Decay rate in steady state = Pausing(b) + Degradation(c)
    # b = rho_ss * L * kp
    b = rho_ss * L * kp
    rate_steady = b + c
    
    # Expected time = P_survive * (1/rate)
    time_steady = prob_survive_filling / rate_steady
    
    load_steady = (rho_ss * L) * time_steady


    # --- 6. Final Density Contribution ---
    # Total Particle-Seconds / (L * tau)
    # Note: We divide by L to get density, and by tau to average over lifetime.
    
    total_load = load_filling + load_steady
    
    rho_0_contrib = total_load / (L * tau)
    
    return rho_0_contrib
end

"""
    calculate_state_densities(L, kp, km, l_ribosome, init, elong, tau)

Decomposes the system into (P0, rho0, Pp, rhoP).
Uses DIRECT Conditional Calculation for rhoP to avoid numerical instability at low kp.
"""
function calculate_state_densities(L, kp, km, l_ribosome, init, elong, tau)
    # --- 1. Calculate Time Fractions (P0, Pp) ---
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    Tp_total = calculate_Tp_total(L, kp, km, l_ribosome, init, elong, tau)
    
    total_life = T0 + Tp_total
    
    if total_life > 1e-20
        P0 = T0 / total_life
        Pp = Tp_total / total_life
    else
        return 0.0, 0.0, 0.0, 0.0
    end

    # --- 2. Unpaused Density (rho0) ---
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        rho0 = 1 / (l_ribosome + sqrt(l_ribosome)) 
    else
        rho0 = init / (elong + init * (l_ribosome - 1))
    end

    # --- 3. Paused Density (rhoP) - DIRECT METHOD ---
    # We calculate Load and Duration separately to avoid (0/0) errors.
    
    if kp <= 1e-20
        rhoP = 0.0
    else
        # A. Calculate the Duration of a SINGLE Jam Event (Tp_dwell)
        # ---------------------------------------------------------
        avg_pos = calculate_avg_pause_position_robust(init, elong, L, kp, tau; ell=l_ribosome)
        N_max = avg_pos / l_ribosome
        
        # Batch sizes
        d = (km / kp); B = d + 1.0 
        Q_max = ceil(Int, N_max / B)
        
        # Effective Step Duration
        dt_eff = tau * (1.0 - exp(-1.0 / (km * tau)))
        p_step = km / (km + 1.0/tau)
        
        # Geometric Sum for Duration
        if abs(1.0 - p_step) < 1e-6
            S_time = Float64(Q_max)
        else
            S_time = (1.0 - p_step^Q_max) / (1.0 - p_step)
        end
        
        # This is the expected duration of the jam GIVEN it exists
        Tp_dwell_conditional = dt_eff * S_time

        # B. Calculate the Total Load of a SINGLE Jam Event (Load_conditional)
        # ------------------------------------------------------------------
        # Re-using the jam logic but without P_entry weighting
        
        # Parameters
        rho_ss = rho0 # Unpaused density
        J_ss   = rho_ss * (elong - init) # Approx flux
        N_start = rho_ss * avg_pos
        
        # Growth Phase Load
        t_grow = 0.0
        if J_ss > 1e-9; t_grow = (N_max - N_start) / J_ss; end
        if t_grow < 0; t_grow = 0.0; end
        
        E_g = exp(-t_grow / tau)
        
        term_linear = tau * (1.0 - E_g) * N_start
        term_growth = J_ss * ((tau^2 * (1.0 - E_g)) - (tau * t_grow * E_g))
        L_growth = term_linear + term_growth
        
        # Dissolution Phase Load
        if abs(1.0 - p_step) < 1e-6
            S1 = Float64(Q_max)
            S2 = Q_max * (Q_max - 1.0) / 2.0
        else
            pQ = p_step^Q_max
            inv_1_p = 1.0 / (1.0 - p_step)
            S1 = (1.0 - pQ) * inv_1_p
            term_a = p_step * (1.0 - p_step^(Q_max - 1)) * (inv_1_p^2)
            term_b = (Q_max - 1) * pQ * inv_1_p
            S2 = term_a - term_b
        end
        L_diss = dt_eff * (N_max * S1 - B * S2)
        
        # Total Load (Particle-Seconds) for this single event
        Total_Load_Conditional = L_growth #+ (E_g * L_diss)
        
        
        # C. Final Calculation
        # rhoP = (Total Particle-Seconds) / (Duration of Event * Length)
        
        if Tp_dwell_conditional > 1e-20
            rhoP = Total_Load_Conditional / (Tp_dwell_conditional * L)
        else
            rhoP = 0.0
        end
    end

    return P0, rho0, Pp, rhoP
end

function calculate_state_densities_corrected(L, kp, km, l_ribosome, init, elong, tau)
    # 1. Weights
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    Tp_total = calculate_Tp_total(L, kp, km, l_ribosome, init, elong, tau)
    
    total_life = T0 + Tp_total
    if total_life <= 1e-20; return 0.0, 0.0, 0.0, 0.0; end
    
    P0 = T0 / total_life
    Pp = Tp_total / total_life
    
    # 2. Conditional rho0 (Unpaused)
    # We use the integral function to get the correct "Empty Start" average density
    # But we must normalize it by T0, not tau.
    
    # Get weighted density contribution (normalized by tau)
    rho0_weighted_contrib = calc_rho_unpaused_contribution(init, l_ribosome, elong, L, kp, tau)
    
    # Convert to conditional density: (Contribution * tau) / T0
    if T0 > 1e-20
        # Check for P0=0 edge case to avoid instability
        rho0 = (rho0_weighted_contrib * tau) / T0
    else
        rho0 = 0.0
    end
    
    # 3. Conditional rhoP (Paused)
    # Uses the direct load/duration logic we fixed previously
    if kp <= 1e-20
        rhoP = 0.0
    else
        # ... [Insert the "DIRECT METHOD" logic for rhoP from the previous response] ...
        # (The logic where we calculate Load_Conditional / Duration_Conditional)
        
        # For brevity, I will call the function wrapper if you have it, 
        # or you can paste the block from the previous "DIRECT METHOD" response.
        # It is crucial to use the version that returns Load / (L * Duration).
        
        # Let's assume you use the logic that worked for you in Method A:
        # Re-calculating Duration and Load directly:
        
        avg_pos = calculate_avg_pause_position_robust(init, elong, L, kp, tau; ell=l_ribosome)
        N_max = avg_pos / l_ribosome
        d = (km / kp); B = d + 1.0; Q_max = ceil(Int, N_max / B)
        dt_eff = tau * (1.0 - exp(-1.0 / (km * tau))); p_step = km / (km + 1.0/tau)
        
        if abs(1.0 - p_step) < 1e-6; S_time = Float64(Q_max)
        else; S_time = (1.0 - p_step^Q_max) / (1.0 - p_step); end
        Tp_dwell = dt_eff * S_time
        
        # Load (Jam + Runoff)
        # We need the TOTAL Paused Load (Jam + Runoff).
        # We can use 'calc_density_paused_total' which returns P_entry * (Load / L*tau)
        
        rhoP_weighted = calc_density_paused_total(km, kp, l_ribosome, L, init, elong, tau)
        
        # Convert to conditional:
        # P_entry is roughly Pp * (tau / Tp_dwell) ... this conversion is tricky.
        # It is safer to calculate Load directly. 
        
        # If 'calc_density_paused_total' is working well for the shape, use it:
        # rhoP = (rhoP_weighted * tau) / Tp_total
        # Note: Use Tp_total (total time paused per life) not Tp_dwell (single event).
        
        if Tp_total > 1e-20
             rhoP = (rhoP_weighted * tau) / Tp_total
        else
             rhoP = 0.0
        end
    end

    return P0, rho0, Pp, rhoP
end

"""
    calc_density_growth_and_runoff(L, kp, l_ribosome, init, elong, tau)

Calculates the weighted density contribution (P_entry * rho_cond) for the PAUSED state,
considering ONLY:
1. The Growth Phase (Linear buildup of ribosomes behind the jam).
2. The Run-off Phase (Downstream ribosomes exiting).
It NEGLECTS the Dissolution Phase (particles leaving the cluster).
"""
function calc_density_growth_and_runoff(L, kp, l_ribosome, init, elong, tau)
    # Safety
    if kp <= 0.0 || tau <= 0.0; return 0.0; end

    # 1. Calculate Entry Probability
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    P_entry = 1.0 - (T0 / tau)
    if P_entry <= 1e-20; return 0.0; end

    # 2. Shared Physics Parameters
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        rho_ss = 1 / (l_ribosome + sqrt(l_ribosome))
    else
        rho_ss = init / (elong + init * (l_ribosome - 1))
    end
    
    v = elong - init; if v <= 1e-9; v = elong; end
    J_ss = rho_ss * v

    # Average Pausing Position
    avg_pos = calculate_avg_pause_position_robust(init, elong, L, kp, tau; ell=l_ribosome)
    
    # ---------------------------------------------------------
    # PART A: GROWTH PHASE LOAD
    # N(t) = N_start + J_ss * t
    # ---------------------------------------------------------
    
    # N_start: Ribosomes already in the segment [0, avg_pos]
    N_start = rho_ss * avg_pos
    
    # N_max: Maximum capacity of the cluster
    N_max = avg_pos / l_ribosome
    
    # Time to fill the cluster (Growth Duration)
    t_grow = 0.0
    if J_ss > 1e-9
        t_grow = (N_max - N_start) / J_ss
    end
    # Clamp: Growth cannot be negative (if already full)
    if t_grow < 0.0; t_grow = 0.0; end
    
    # Integration: Int_0^t_grow (N_start + J*t) * exp(-t/tau) dt
    
    E_g = exp(-t_grow / tau)
    
    # 1. Constant Term: Int N_start * exp(-t/tau)
    # = N_start * tau * (1 - E_g)
    load_const = N_start * tau * (1.0 - E_g)
    
    # 2. Linear Term: Int J * t * exp(-t/tau)
    # Int t*e^-ct = (1/c^2)(1 - e^-cT) - (T/c)e^-cT
    # Here c = 1/tau, so 1/c = tau, 1/c^2 = tau^2
    term_linear_integral = (tau^2 * (1.0 - E_g)) - (tau * t_grow * E_g)
    load_linear = J_ss * term_linear_integral
    
    Total_Load_Growth = load_const + load_linear


    # ---------------------------------------------------------
    # PART B: RUN-OFF LOAD
    # Downstream particles exiting (Same as before)
    # ---------------------------------------------------------
    
    wall_pos = 2.0 * avg_pos
    
    if wall_pos < L
        stream_len = avg_pos
        gap_len = L - wall_pos
    else
        stream_len = L - avg_pos
        gap_len = 0.0
    end
    if stream_len < 0; stream_len = 0.0; end

    t_start_exit = gap_len / v
    t_end_exit   = (gap_len + stream_len) / v
    
    # Integral 1: Gap Delay (Constant Mass)
    load_gap = 0.0
    if t_start_exit > 1e-9
        mass_gap = rho_ss * stream_len
        integral_gap = tau * (1.0 - exp(-t_start_exit / tau))
        load_gap = mass_gap * integral_gap
    end
    
    # Integral 2: Exit Phase (Linear Decay)
    load_exit = 0.0
    if t_end_exit > t_start_exit
        duration = t_end_exit - t_start_exit
        prob_survive_gap = exp(-t_start_exit / tau)
        
        # Constant Part
        int_const_exit = tau * (1.0 - exp(-duration / tau))
        term_1_exit = (rho_ss * stream_len) * int_const_exit
        
        # Decay Part (Mass leaving)
        # Int (v * u) * exp(-u/tau)
        val_upper = exp(-duration/tau) * (-duration/tau - 1.0)
        val_lower = -1.0
        term_2_exit = (rho_ss * v) * (tau^2 * (val_upper - val_lower))
        
        load_exit = prob_survive_gap * (term_1_exit - term_2_exit)
    end
    
    Total_Load_Runoff = load_gap + load_exit

    # ---------------------------------------------------------
    # COMBINE
    # ---------------------------------------------------------
    
    # Total Load conditional on entering paused state
    Total_Load_Conditional = Total_Load_Growth + Total_Load_Runoff
    
    # Weighted Density Contribution = P_entry * (Load / (L * tau))
    # Note: We divide by tau here to get the time-averaged density contribution over the lifetime.
    
    rho_weighted_contribution = P_entry * (Total_Load_Conditional / (L * tau))
    
    return rho_weighted_contribution
end

function calculate_state_densities_run_off_w_growth(L, kp, km, l_ribosome, init, elong, tau)
    # 1. Weights
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    Tp_total = calculate_Tp_total(L, kp, km, l_ribosome, init, elong, tau)
    
    total_life = T0 + Tp_total
    if total_life <= 1e-20; return 0.0, 0.0, 0.0, 0.0; end
    
    P0 = T0 / total_life
    Pp = Tp_total / total_life
    
    # 2. Conditional rho0 (Unpaused)
    rho0_weighted_contrib = calc_rho_unpaused_contribution(init, l_ribosome, elong, L, kp, tau)
    
    if T0 > 1e-20
        rho0 = (rho0_weighted_contrib * tau) / T0
    else
        rho0 = 0.0
    end
    
    # 3. Conditional rhoP (Paused) - GROWTH + RUNOFF ONLY
    if kp <= 1e-20
        rhoP = 0.0
    else
        # Calculate the weighted contribution (Growth + Runoff)
        rhoP_weighted = calc_density_growth_and_runoff(L, kp, l_ribosome, init, elong, tau)
        
        # Normalize to get conditional density
        # rhoP = (Contribution * tau) / Tp_total
        if Tp_total > 1e-20
             rhoP = (rhoP_weighted * tau) / Tp_total
        else
             rhoP = 0.0
        end
    end

    return P0, rho0, Pp, rhoP
end

"""
    calculate_state_densities_corrected(L, kp, km, l_ribosome, init, elong, tau)

Calculates (P0, rho0, Pp, rhoP).
Uses a CONTINUOUS INTEGRATION approach for the Paused Density (rhoP), considering:
1. Growth Phase (Ramp up)
2. Dissolution Phase (Jam decay) -> Crucial for short lattices!
3. Run-off Phase (Exit stream)
"""
function calculate_state_densities_corrected2(L, kp, km, l_ribosome, init, elong, tau)
    # --- 1. Weights & Probabilities ---
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    Tp_total = calculate_Tp_total(L, kp, km, l_ribosome, init, elong, tau)
    
    total_life = T0 + Tp_total
    if total_life <= 1e-20; return 0.0, 0.0, 0.0, 0.0; end
    
    P0 = T0 / total_life
    Pp = Tp_total / total_life
    
    # Entry Probability (used to weight the paused load)
    P_entry = 1.0 - (T0 / tau)
    if P_entry < 0; P_entry = 0.0; end


    # --- 2. Conditional rho0 (Unpaused) ---
    # Get weighted density contribution (normalized by tau)
    rho0_weighted_contrib = calc_rho_unpaused_contribution(init, l_ribosome, elong, L, kp, tau)
    
    if T0 > 1e-20
        # Convert to conditional density: (Contribution * tau) / T0
        rho0 = (rho0_weighted_contrib * tau) / T0
    else
        rho0 = 0.0
    end


    # --- 3. Conditional rhoP (Paused) ---
    # We calculate the Total Load (Growth + Dissolution + Runoff) 
    # and normalize by the Total Paused Time (Tp_total).
    
    if kp <= 1e-20 || P_entry <= 1e-20
        rhoP = 0.0
    else
        # --- Shared Physics Parameters ---
        if init/elong >= 1/ (1 + sqrt(l_ribosome))
            rho_ss = 1 / (l_ribosome + sqrt(l_ribosome))
        else
            rho_ss = init / (elong + init * (l_ribosome - 1))
        end
        v = elong - init; if v <= 1e-9; v = elong; end
        J_ss = rho_ss * v
        
        avg_pos = calculate_avg_pause_position_robust(init, elong, L, kp, tau; ell=l_ribosome)
        
        # --- PART A: GROWTH PHASE (Ramp Up) ---
        N_start = rho_ss * avg_pos
        N_max = avg_pos / l_ribosome
        
        t_grow = 0.0
        if J_ss > 1e-9; t_grow = (N_max - N_start) / J_ss; end
        if t_grow < 0; t_grow = 0.0; end
        
        E_g = exp(-t_grow / tau)
        
        # Load = Int_0^t_grow (N_start + J*t)*e^-t/tau
        load_growth = N_start * tau * (1.0 - E_g) + 
                      J_ss * ((tau^2 * (1.0 - E_g)) - (tau * t_grow * E_g))
        
        # --- PART B: DISSOLUTION PHASE (The Jam) ---
        # Continuous approximation: Linear decay from N_max to 0 over T_clear
        B = (km / kp) + 1.0
        T_clear_theoretical = (N_max / B) / km
        
        if T_clear_theoretical > 1e-9
            E_d = exp(-T_clear_theoretical / tau)
            
            # Int_0^T (N_max - N_max/T * u) * e^-u/tau
            term_const = N_max * tau * (1.0 - E_d)
            
            val_upper = E_d * (-T_clear_theoretical/tau - 1.0)
            val_lower = -1.0
            term_linear = (N_max / T_clear_theoretical) * (tau^2 * (val_upper - val_lower))
            
            load_dissolution_integral = term_const - term_linear
            
            # Weighted by survival of growth phase
            load_dissolution = E_g * load_dissolution_integral
        else
            load_dissolution = 0.0
        end

        # --- PART C: RUN-OFF PHASE (Exit Stream) ---
        wall_pos = 2.0 * avg_pos
        if wall_pos < L
            stream_len = avg_pos; gap_len = L - wall_pos
        else
            stream_len = L - avg_pos; gap_len = 0.0
        end
        if stream_len < 0; stream_len = 0.0; end

        t_start_exit = gap_len / v
        t_end_exit   = (gap_len + stream_len) / v
        
        # Gap Delay Load
        load_gap = 0.0
        if t_start_exit > 1e-9
            mass_gap = rho_ss * stream_len
            integral_gap = tau * (1.0 - exp(-t_start_exit / tau))
            load_gap = mass_gap * integral_gap
        end
        
        # Exit Decay Load
        load_exit = 0.0
        if t_end_exit > t_start_exit
            duration = t_end_exit - t_start_exit
            prob_survive_gap = exp(-t_start_exit / tau)
            
            int_const_exit = tau * (1.0 - exp(-duration / tau))
            term_1_exit = (rho_ss * stream_len) * int_const_exit
            
            val_upper = exp(-duration/tau) * (-duration/tau - 1.0)
            val_lower = -1.0
            term_2_exit = (rho_ss * v) * (tau^2 * (val_upper - val_lower))
            
            load_exit = prob_survive_gap * (term_1_exit - term_2_exit)
        end
        
        Total_Load_Runoff = load_gap + load_exit

        # --- FINAL CALCULATION ---
        Total_Load_Conditional = load_growth  + Total_Load_Runoff + load_dissolution
        
        # Weighted Density Contribution (normalized by tau)
        rhoP_weighted_contrib = P_entry * (Total_Load_Conditional / (L * tau))
        
        # Convert to conditional density: (Contribution * tau) / Tp_total
        # This effectively calculates: Total_Load_Conditional / (L * Tp_total_effective)
        if Tp_total > 1e-20
             rhoP = (rhoP_weighted_contrib * tau) / Tp_total
        else
             rhoP = 0.0
        end
    end

    return P0, rho0, Pp, rhoP
end

function calculate_state_densities_corrected_plateau(L, kp, km, l_ribosome, init, elong, tau)
    # --- 1. Weights & Probabilities ---
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    Tp_total = calculate_Tp_total(L, kp, km, l_ribosome, init, elong, tau)
    
    total_life = T0 + Tp_total
    if total_life <= 1e-20; return 0.0, 0.0, 0.0, 0.0; end
    
    P0 = T0 / total_life
    Pp = Tp_total / total_life
    
    P_entry = 1.0 - (T0 / tau)
    if P_entry < 0; P_entry = 0.0; end

    # --- 2. Conditional rho0 (Unpaused) ---
    rho0_weighted_contrib = calc_rho_unpaused_contribution(init, l_ribosome, elong, L, kp, tau)
    rho0 = (T0 > 1e-20) ? (rho0_weighted_contrib * tau) / T0 : 0.0

    # --- 3. Conditional rhoP (Paused) ---
    if kp <= 1e-20 || P_entry <= 1e-20
        rhoP = 0.0
    else
        # --- Physics Parameters ---
        if init/elong >= 1/ (1 + sqrt(l_ribosome))
            rho_ss = 1 / (l_ribosome + sqrt(l_ribosome))
        else
            rho_ss = init / (elong + init * (l_ribosome - 1))
        end
        v = elong - init; if v <= 1e-9; v = elong; end
        J_ss = rho_ss * v
        avg_pos = calculate_avg_pause_position_robust(init, elong, L, kp, tau; ell=l_ribosome)
        
        # --- PART A: GROWTH PHASE (Ramp Up) ---
        N_start = rho_ss * avg_pos
        N_max = avg_pos / l_ribosome
        
        t_grow = 0.0
        if J_ss > 1e-9; t_grow = (N_max - N_start) / J_ss; end
        if t_grow < 0; t_grow = 0.0; end
        E_g = exp(-t_grow / tau)
        
        load_growth = N_start * tau * (1.0 - E_g) + 
                      J_ss * ((tau^2 * (1.0 - E_g)) - (tau * t_grow * E_g))
        
        # --- PART B: PLATEAU PHASE (Static Jam) ---
        # The jam sits at N_max for the time it takes the antibiotic to unbind.
        # We neglect the subsequent shrinking. 
        # T_wait = 1/km
        
        T_wait_theoretical = 1.0 / km
        
        # Effective duration limited by degradation
        t_wait_eff = tau * (1.0 - exp(-T_wait_theoretical / tau))
        
        # Load = Mass * Duration
        # Weighted by the probability that the system survived growth (E_g)
        load_plateau = E_g * (N_max * t_wait_eff)

        # --- PART C: RUN-OFF PHASE (Exit Stream) ---
        wall_pos = 2.0 * avg_pos
        if wall_pos < L; stream_len = avg_pos; gap_len = L - wall_pos
        else; stream_len = L - avg_pos; gap_len = 0.0; end
        if stream_len < 0; stream_len = 0.0; end

        t_start_exit = gap_len / v; t_end_exit = (gap_len + stream_len) / v
        
        # Gap Delay Load
        load_gap = 0.0
        if t_start_exit > 1e-9
            mass_gap = rho_ss * stream_len
            integral_gap = tau * (1.0 - exp(-t_start_exit / tau))
            load_gap = mass_gap * integral_gap
        end
        
        # Exit Decay Load
        load_exit = 0.0
        if t_end_exit > t_start_exit
            duration = t_end_exit - t_start_exit
            prob_survive_gap = exp(-t_start_exit / tau)
            int_const_exit = tau * (1.0 - exp(-duration / tau))
            term_1_exit = (rho_ss * stream_len) * int_const_exit
            val_upper = exp(-duration/tau) * (-duration/tau - 1.0); val_lower = -1.0
            term_2_exit = (rho_ss * v) * (tau^2 * (val_upper - val_lower))
            load_exit = prob_survive_gap * (term_1_exit - term_2_exit)
        end
        Total_Load_Runoff = load_gap + load_exit

        # --- FINAL CALCULATION ---
        # Summing Growth + Plateau + Runoff
        Total_Load_Conditional = load_growth + load_plateau + Total_Load_Runoff
        
        rhoP_weighted_contrib = P_entry * (Total_Load_Conditional / (L * tau))
        
        if Tp_total > 1e-20
             rhoP = (rhoP_weighted_contrib * tau) / Tp_total
        else
             rhoP = 0.0
        end
    end

    return P0, rho0, Pp, rhoP
end

function calculate_state_densities_corrected_platea2(L, kp, km, l_ribosome, init, elong, tau)
    # --- 1. Weights & Probabilities ---
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    Tp_total = calculate_Tp_total(L, kp, km, l_ribosome, init, elong, tau)
    
    total_life = T0 + Tp_total
    if total_life <= 1e-20; return 0.0, 0.0, 0.0, 0.0; end
    
    P0 = T0 / total_life
    Pp = Tp_total / total_life
    
    P_entry = 1.0 - (T0 / tau)
    if P_entry < 0; P_entry = 0.0; end

    # --- 2. Conditional rho0 (Unpaused) ---
    rho0_weighted_contrib = calc_rho_unpaused_contribution(init, l_ribosome, elong, L, kp, tau)
    rho0 = (T0 > 1e-20) ? (rho0_weighted_contrib * tau) / T0 : 0.0

    # --- 3. Conditional rhoP (Paused) ---
    if kp <= 1e-20 || P_entry <= 1e-20
        rhoP = 0.0
    else
        # --- Physics Parameters ---
        if init/elong >= 1/ (1 + sqrt(l_ribosome))
            rho_ss = 1 / (l_ribosome + sqrt(l_ribosome))
        else
            rho_ss = init / (elong + init * (l_ribosome - 1))
        end
        v = elong - init; if v <= 1e-9; v = elong; end
        J_ss = rho_ss * v
        avg_pos = calculate_avg_pause_position_robust(init, elong, L, kp, tau; ell=l_ribosome)
        
        # === FIX 1: DISCRETE JAM CAPACITY ===
        # Continuous: N_max = 1.8 -> Overestimation.
        # Discrete:   N_max = floor(1.8) = 1.0 -> Correct.
        
        N_max_continuous = avg_pos / l_ribosome
        N_max = floor(N_max_continuous)
        
        # N_start must also be capped by the physical capacity
        N_start = min(N_max, rho_ss * avg_pos)
        
        # --- PART A: GROWTH PHASE (Ramp Up) ---
        t_grow = 0.0
        # Only grow if we aren't already full
        if J_ss > 1e-9 && N_max > N_start
             t_grow = (N_max - N_start) / J_ss
        end
        if t_grow < 0; t_grow = 0.0; end
        
        E_g = exp(-t_grow / tau)
        
        load_growth = N_start * tau * (1.0 - E_g) + 
                      J_ss * ((tau^2 * (1.0 - E_g)) - (tau * t_grow * E_g))
        
        # --- PART B: PLATEAU PHASE (Static Jam) ---
        # The jam sits at N_max until the antibiotic unbinds (or mRNA dies).
        T_wait_theoretical = 1.0 / km
        
        # === FIX 2: TIME CORRECTION ===
        # The Plateau only starts AFTER the growth is finished.
        T_plateau_duration = max(0.0, T_wait_theoretical - t_grow)
        
        # Effective duration limited by degradation
        t_wait_eff = tau * (1.0 - exp(-T_plateau_duration / tau))
        
        # Load = Mass * Duration * Probability of surviving growth (E_g)
        load_plateau = E_g * (N_max * t_wait_eff)

        # --- PART C: RUN-OFF PHASE (Exit Stream) ---
        wall_pos = 2.0 * avg_pos
        
        # === FIX 3: SPACE CORRECTION ===
        # Subtract the paused particle from the stream to avoid double-counting.
        particle_len = (rho_ss > 1e-9) ? (1.0 / rho_ss) : 0.0
        
        if wall_pos < L
            stream_len = max(0.0, avg_pos - particle_len)
            gap_len = L - wall_pos
        else
            stream_len = max(0.0, (L - avg_pos) - particle_len)
            gap_len = 0.0
        end

        t_start_exit = gap_len / v
        t_end_exit = (gap_len + stream_len) / v
        
        # Gap Delay Load
        load_gap = 0.0
        if t_start_exit > 1e-9
            mass_gap = rho_ss * stream_len
            integral_gap = tau * (1.0 - exp(-t_start_exit / tau))
            load_gap = mass_gap * integral_gap
        end
        
        # Exit Decay Load
        load_exit = 0.0
        if t_end_exit > t_start_exit
            duration = t_end_exit - t_start_exit
            prob_survive_gap = exp(-t_start_exit / tau)
            int_const_exit = tau * (1.0 - exp(-duration / tau))
            term_1_exit = (rho_ss * stream_len) * int_const_exit
            val_upper = exp(-duration/tau) * (-duration/tau - 1.0); val_lower = -1.0
            term_2_exit = (rho_ss * v) * (tau^2 * (val_upper - val_lower))
            load_exit = prob_survive_gap * (term_1_exit - term_2_exit)
        end
        Total_Load_Runoff = load_gap + load_exit

        # --- FINAL CALCULATION ---
        Total_Load_Conditional = load_growth + load_plateau + Total_Load_Runoff
        
        rhoP_weighted_contrib = P_entry * (Total_Load_Conditional / (L * tau))
        
        if Tp_total > 1e-20
             rhoP = (rhoP_weighted_contrib * tau) / Tp_total
        else
             rhoP = 0.0
        end
    end

    return P0, rho0, Pp, rhoP
end

function calculate_state_densities_corrected_plateau4(L, kp, km, l_ribosome, init, elong, tau)
    # --- 1. Weights & Probabilities ---
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    Tp_total = calculate_Tp_total(L, kp, km, l_ribosome, init, elong, tau)
    
    total_life = T0 + Tp_total
    if total_life <= 1e-20; return 0.0, 0.0, 0.0, 0.0; end
    
    P0 = T0 / total_life
    Pp = Tp_total / total_life
    P_entry = max(0.0, 1.0 - (T0 / tau))

    # --- 2. Conditional rho0 (Unpaused) ---
    rho0_weighted_contrib = calc_rho_unpaused_contribution(init, l_ribosome, elong, L, kp, tau)
    rho0 = (T0 > 1e-20) ? (rho0_weighted_contrib * tau) / T0 : 0.0

    # --- 3. Conditional rhoP (Paused) ---
    if kp <= 1e-20 || P_entry <= 1e-20
        rhoP = 0.0
    else
        # ... [Physics Parameters] ...
        if init/elong >= 1/ (1 + sqrt(l_ribosome))
            rho_ss = 1 / (l_ribosome + sqrt(l_ribosome))
        else
            rho_ss = init / (elong + init * (l_ribosome - 1))
        end
        v = elong - init; if v <= 1e-9; v = elong; end
        J_ss = rho_ss * v
        avg_pos = calculate_avg_pause_position_robust(init, elong, L, kp, tau; ell=l_ribosome)
        
        # Discrete Jam Capacity
        N_max = floor(avg_pos / l_ribosome)
        N_start = min(N_max, rho_ss * avg_pos)
        
        # Combined Decay Rate (Death + Unbinding)
        # This is the rate at which the "Growing/Full Jam" state ends.
        k_eff = (1.0 / tau) + km

        # --- PART A: GROWTH PHASE ---
        t_grow = 0.0
        if J_ss > 1e-9 && N_max > N_start
             t_grow = (N_max - N_start) / J_ss
        end
        if t_grow < 0; t_grow = 0.0; end
        
        # Probability of surviving growth without dying OR unbinding
        E_g_comb = exp(-t_grow * k_eff)
        
        # Growth Load (Integral against k_eff)
        term_const = N_start * (1.0 - E_g_comb) / k_eff
        term_linear_int = ((1.0 / k_eff^2) * (1.0 - E_g_comb)) - ((t_grow / k_eff) * E_g_comb)
        load_growth = term_const + (J_ss * term_linear_int)
        
        # --- PART B: PLATEAU PHASE ---
        # The jam sits at N_max, decaying at rate k_eff.
        # It starts at t_grow (prob E_g_comb) and lasts until Death or Unbinding.
        # Integral_0^inf (N_max * exp(-k_eff * t)) = N_max / k_eff
        
        load_plateau = E_g_comb * (N_max / k_eff)

        # --- PART C: DISSOLUTION PHASE (The Correction) ---
        # If unbinding happens (rate km), the cluster shrinks.
        # Prob of entering Dissolution = P(Survive Growth) * P(Unbind vs Death)
        prob_unbind = km / k_eff
        weight_dissolution = E_g_comb * prob_unbind
        
        if weight_dissolution > 1e-9
             # Standard Discrete Dissolution Logic (Shrinking Batches)
             d = (km / kp); B = d + 1.0; Q_max = ceil(Int, N_max / B)
             
             # Effective step duration (limited by mRNA lifetime tau)
             # Note: Once unbound, only tau matters (km is done).
             dt_step = tau * (1.0 - exp(-1.0 / (km * tau)))
             p_step = km / (km + 1.0/tau)
             
             if abs(1.0 - p_step) < 1e-6
                 S1 = Float64(Q_max)
                 S2 = Q_max * (Q_max - 1.0) / 2.0
             else
                 pQ = p_step^Q_max
                 inv_1_p = 1.0 / (1.0 - p_step)
                 S1 = (1.0 - pQ) * inv_1_p
                 term_a = p_step * (1.0 - p_step^(Q_max - 1)) * (inv_1_p^2)
                 term_b = (Q_max - 1) * pQ * inv_1_p
                 S2 = term_a - term_b
             end
             
             # Total Load of the shrinking process
             load_dissolution_raw = dt_step * (N_max * S1 - B * S2)
             
             load_dissolution = weight_dissolution * load_dissolution_raw
        else
             load_dissolution = 0.0
        end

        # --- PART D: RUN-OFF PHASE ---
        # (Standard logic with Space Correction)
        wall_pos = 2.0 * avg_pos
        particle_len = (rho_ss > 1e-9) ? (1.0 / rho_ss) : 0.0
        if wall_pos < L; stream_len = max(0.0, avg_pos - particle_len); gap_len = L - wall_pos
        else; stream_len = max(0.0, (L - avg_pos) - particle_len); gap_len = 0.0; end

        t_start_exit = gap_len / v; t_end_exit = (gap_len + stream_len) / v
        
        load_gap = 0.0
        if t_start_exit > 1e-9
            mass_gap = rho_ss * stream_len
            integral_gap = tau * (1.0 - exp(-t_start_exit / tau))
            load_gap = mass_gap * integral_gap
        end
        load_exit = 0.0
        if t_end_exit > t_start_exit
            duration = t_end_exit - t_start_exit
            prob_survive_gap = exp(-t_start_exit / tau)
            int_const_exit = tau * (1.0 - exp(-duration / tau))
            term_1_exit = (rho_ss * stream_len) * int_const_exit
            val_upper = exp(-duration/tau) * (-duration/tau - 1.0); val_lower = -1.0
            term_2_exit = (rho_ss * v) * (tau^2 * (val_upper - val_lower))
            load_exit = prob_survive_gap * (term_1_exit - term_2_exit)
        end
        Total_Load_Runoff = load_gap + load_exit

        # --- FINAL ---
        Total_Load_Conditional = load_growth + load_plateau + load_dissolution + Total_Load_Runoff
        rhoP_weighted_contrib = P_entry * (Total_Load_Conditional / (L * tau))
        
        if Tp_total > 1e-20
             rhoP = (rhoP_weighted_contrib * tau) / Tp_total
        else
             rhoP = 0.0
        end
    end
    return P0, rho0, Pp, rhoP
end


using SpecialFunctions

# --- 1. Helper Functions ---


"""
    calculate_rho0_appendix(L, kp, tau, T0)

Calculates the time-averaged density of the unpaused state (rho_0)
using the exact analytical result derived in Appendix A.

Formula: rho0 = (1 / L*kp) * ( (1/T0) - (1/tau) )

# Arguments
- `L`: Lattice length
- `kp`: Pausing rate
- `tau`: mRNA lifetime
- `T0`: Expected duration of the unpaused state (calculated via error functions)
"""
function calculate_rho0_appendix(L, kp, tau, l_ribosome, init, elong)
    # --- Stability Check ---
    # As kp -> 0, T0 -> tau. The term becomes (1/tau - 1/tau) / 0 = 0/0.
    # Physically, this represents a standard TASEP without pausing.
    # To prevent NaN, we check for small kp.
    if kp <= 1e-12
        # In this limit, rho0 is simply the time-averaged density of a 
        # filling TASEP under degradation, which requires the explicit integral.
        # Alternatively, return 0.0 or the steady-state rho_ss if T0 ~ tau.
        return 0.0 
    end
    T0 = calculate_T0(init, elong, L, kp, tau; ell=l_ribosome)
    # --- Analytical Formula (Appendix A) ---
    # Term A: 1 / T0 (Effective rate of leaving unpaused state)
    rate_exit_total = 1.0 / T0
    
    # Term B: 1 / tau (Rate of leaving due to degradation)
    rate_exit_decay = 1.0 / tau
    
    # Difference normalized by total pausing capacity (L * kp)
    rho0 = (1.0 / (L * kp)) * (rate_exit_total - rate_exit_decay)
    
    return rho0
end


using SpecialFunctions

function get_tasep_params(alpha, epsilon, ell)
    if alpha/epsilon >= 1/ (1 + sqrt(ell))
        rho = 1 / (ell + sqrt(ell))
        J = epsilon / (1 + sqrt(ell))^2
    else
        rho = alpha / (epsilon + alpha * (ell - 1))
        J = (alpha * (epsilon - alpha)) / (epsilon + alpha * (ell - 1))
    end
    v = epsilon - alpha
    return rho, J, v
end

"""
    calculate_rho_p_unified(L, kp, km, ell, alpha, epsilon, tau; include_gap_travel=true)

Calculates the TOTAL Paused Density by summing four physical phases:
1. Jam Formation & Persistence (Decays at tau_eff)
2. Run-off (Downstream particles exiting, Decays at tau)
3. Uncoiling (The cluster melting after unbinding)
4. Recovery (The lattice refilling after melting)
"""
function calculate_rho_p_unified(L, kp, km, ell, alpha, epsilon, tau; include_gap_travel=true)
    rho_ss, J_ss, v = get_tasep_params(alpha, epsilon, ell)
    
    # 1. Timescales
    # tau_eff: Average time until Jam ends (Death OR Unbinding)
    k_eff = (1.0 / tau) + km
    #tau_eff = 1.0 / k_eff
    tau_eff = tau
    # 2. Geometry
    if kp <= 1e-12; Xf = L/2.0; else
        Xf = calculate_avg_pause_position_robust(alpha, epsilon, L, kp, tau; ell=ell)
    end
    N_init = rho_ss * Xf
    N_max = Xf / ell
    
    # --- PHASE 1: THE JAM (Static) ---
    # Calculates density while the antibiotic is bound.
    t_cf = (J_ss > 1e-9 && N_max > N_init) ? (N_max - N_init)/J_ss : 0.0
    E_cf = exp(-t_cf / tau_eff) # Prob of surviving formation
    
    # Load during formation + persistence
    L_cf = (N_init * tau_eff * (1 - E_cf)) + 
           (J_ss * (tau_eff^2 * (1 - E_cf) - tau_eff * t_cf * E_cf))
    L_c  = N_max * tau_eff * E_cf
    
    L_Jam = L_cf + L_c

    # --- PHASE 2: RUN-OFF (Downstream) ---
    # Calculates density of particles downstream of the jam exiting.
    # Decays with 'tau' (mRNA death), unaffected by unbinding.
    Xs = 2.0 * Xf
    L_ro = 0.0
    if Xs >= L 
        t_ro = (L - Xf) / v
        E_ro = exp(-t_ro / tau)
        L_ro = (rho_ss * (L - Xf) * tau) - (J_ss * tau^2 * (1.0 - E_ro))
    else
        t_start = (L - Xs) / v; t_exit = (L - Xf) / v
        E_start = exp(-t_start/tau); E_ro = exp(-(t_exit-t_start)/tau)
        
        if include_gap_travel
            L_ro += (rho_ss * Xf) * tau * (1.0 - E_start)
        end
        L_ro += E_start * ((rho_ss * Xf * tau) - (J_ss * tau^2 * (1.0 - E_ro)))
    end

    # --- PHASE 3: UNCOILING (Dissolution) ---
    # If unbinding happens, the cluster melts.
    # Prob of Unbinding = km / k_eff
    p_unbind = km / k_eff
    
    # Time to clear N particles: N * (ell/v)
    # Load = Triangle area = 0.5 * N_max * Time
    t_melt = N_max * (ell / v)
    load_melt_raw = 0.5 * N_max * t_melt
    
    # Weighted by: Survival of formation (E_cf) * Prob of Unbinding
    L_uncoil = (E_cf * p_unbind) * load_melt_raw

    # --- PHASE 4: RECOVERY (Refill) ---
    # After melting, the lattice refills to rho_ss for the remaining time.
    # Time recovered = tau - tau_eff - t_melt
    # (We subtract t_melt because density is accounted for in Phase 3)
    time_remaining = max(0.0, tau - tau_eff - t_melt)
    
    L_recovery = (E_cf * p_unbind) * (rho_ss * L) * time_remaining

    # --- FINAL SUM ---
    L_total = L_Jam  + L_ro #+ L_uncoil + L_recovery #+ 
    
    return L_total / (L * tau)
end
using QuadGK

"""
    calculate_rho_p_conditional(L, kp, km, ell, alpha, epsilon, tau)

Calculates the Conditional Paused Density (rho_p) defined in the 'Potential Fix' document.
This value represents the density of the ribosome cluster *given that a pause is active*.
It does NOT include the P_p weighting (you must multiply by P_p externally).
"""
using QuadGK

function calculate_rho_p_conditional(L, kp, km, ell, alpha, epsilon, tau)
    # --- 0. Parameters ---
    rho_ss, J_ss, v = get_tasep_params(alpha, epsilon, ell)
    t_L = L / v 
    
    k_eff = (1.0/tau) + km
    tau_eff = tau#1.0 / k_eff

    # --- 1. Probability Functions [cite: 8, 11, 13] ---
    function lambda(t)
        return (t < t_L) ? (J_ss * kp * t) : (rho_ss * L * kp)
    end

    Lambda_at_tL = 0.5 * J_ss * kp * t_L^2
    function Lambda(t)
        return (t < t_L) ? (0.5 * J_ss * kp * t^2) : (Lambda_at_tL + rho_ss * L * kp * (t - t_L))
    end

    F(t) = lambda(t) * exp(-Lambda(t) - (t / tau))

    # --- 2. Load Functions [cite: 24, 33, 36] ---
    
    # HELPER: Calculate correct Jam Formation Time 
    function get_t_cf(Xf)
        if J_ss < 1e-12; return 0.0; end
        # Time to fill = (Capacity - Initial_Load) / Flux
        # Capacity = Xf / ell
        # Initial_Load = rho_ss * Xf
        term = (1.0/ell) - rho_ss
        return (Xf * term) / J_ss
    end

    # A. Cluster Formation Load (During growth) [cite: 24]
    function get_L_cf(Xf)
        t_cf = get_t_cf(Xf)
        
        # Standard integral of linear growth weighted by survival
        E_cf = exp(-t_cf / tau_eff)
        term1 = rho_ss * Xf * tau_eff * (1.0 - E_cf)
        term2 = J_ss * (tau_eff^2 * (1.0 - E_cf) - tau_eff * t_cf * E_cf)
        return term1 + term2
    end

    # B. Cluster Persistence Load (After growth) [cite: 26, 38]
    # THIS WAS MISSING. It accounts for the jam sitting there after it is full.
    function get_L_c(Xf)
        t_cf = get_t_cf(Xf)
        
        # Max capacity of the jam
        N_max = Xf / ell
        
        # Probability the jam actually finishes forming before death/unbinding
        prob_survival = exp(-t_cf / tau_eff)
        
        # If it survives formation, it persists for the remaining lifetime (avg tau_eff)
        # Load = Mass * Time
        return N_max * tau_eff * prob_survival
    end

    # C. Run-off Load (Steady State) [cite: 33]
    function get_L_ro_ss(Xf)
        t_ro = (L - Xf) / v
        return (rho_ss * (L - Xf) * tau) - (J_ss * tau^2 * (1.0 - exp(-t_ro / tau)))
    end

    # D. Run-off Load (Transient) [cite: 36]
    function get_L_ro_transient(Xf)
        if 2.0 * Xf >= L; return get_L_ro_ss(Xf); end
        t_start = (L - 2.0*Xf) / v
        t_ro = Xf / v
        return exp(-t_start / tau) * ((rho_ss * Xf * tau) - (J_ss * tau^2 * (1.0 - exp(-t_ro / tau))))
    end

    # --- 3. Integrate Total Paused Load  ---
    # We sum all three components: Formation + Persistence + Runoff
    
    # Integral A: Filling Phase
    integrand_filling(t) = F(t) * (get_L_cf(v*t/2.0) + get_L_c(v*t/2.0) + get_L_ro_transient(v*t/2.0))
    load_filling, _ = quadgk(integrand_filling, 0.0, t_L)

    # Integral B: Steady State Phase
    # Spatially average the sum of loads
    avg_load_ss, _ = quadgk(x -> get_L_cf(x) + get_L_c(x) + get_L_ro_ss(x), 0.0, L)
    avg_load_ss /= L
    
    prob_ss, _ = quadgk(F, t_L, Inf)
    load_ss = prob_ss * avg_load_ss

    total_integrated_load = load_filling + load_ss

    # --- 4. Normalization by Tp [cite: 48] ---
    a_param = 0.5 * J_ss * kp
    T_fill, _ = quadgk(t -> exp(-a_param * t^2 - t/tau), 0.0, t_L)
    
    num = exp(-a_param * t_L^2 - t_L/tau)
    den = (rho_ss * L * kp) + (1.0/tau)
    T_ss = num / den
    
    T0 = T_fill + T_ss
    Tp = tau - T0

    if Tp < 1e-12; return 0.0; end

    # Return density conditional on being paused
    return total_integrated_load / (L * Tp)
end