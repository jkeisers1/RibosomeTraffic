using Random
using Distributions
using StatsBase
using JLD
#using Plots
print("hello")
include("Basic_functions_obc.jl")
#J0_l(ρ) = (elong*ρ*(1-ρ*l_ribosome)) / (1 + ρ*(l_ribosome-1))
function get_elongation_vector_obc(lattice, rates)
    elongation_vector = zeros(length(lattice))

    for i in 1:length(lattice)-1
        if lattice[i] == 1 && lattice[i+1] == 0
            elongation_vector[i] += rates[i] # adding rate allows for weighted sampling
        end
    end
    
    if lattice[end] == 1  #termination
        elongation_vector[end] += rates[end]
    end
    
    if lattice[2] == 0 # initiation rdy
        elongation_vector[1] += rates[1]
    end

    return elongation_vector
end

function get_internal_state_vec_obc(lattice, k₊, mobile) # where τ = waiting time paused state, f = rate entering mobile state

    internal_state_vec = zeros(length(lattice))

    for i in eachindex(lattice[1:end])

        if lattice[i] == 1
            internal_state_vec[i] = k₊
            mobile[1] += 1
        end
    end
    return internal_state_vec
end

function elongation_process_obc(
    elongation_vector, w_elong, internal_state_vec, 
    posR, lattice, rates, J, k₊, mobile,
    l_ribosome, track_site, jammed, tot_part, num_clust,
    J_c, J_w, paused,c_end_pos_tracker
    )
    #=
    during the elongation process, the propensities have to be updated each time the system evolves.
    Three particles classes:
        paused (with rate 1/τ or k₋)
        jammed (particles jammed by next one)
        mobile (particles able to contribute to the current)
    =#

    moving_particle = sample(posR, w_elong) # which particles moves, the position of the particle is weighted by the rates at the site
    
    
    if moving_particle == c_end_pos_tracker[1]
        c_end_pos_tracker[1] +=1
    end


    if 2 <= moving_particle <= length(lattice)-1#-l_ribosome # bulk hopping

        change_elong_vector_bulk_obc!(
            elongation_vector, moving_particle, lattice, 
            rates, internal_state_vec, k₊, l_ribosome, track_site, 
            jammed, mobile, num_clust
            )
        # changing lattice
        lattice[moving_particle] = 0
        lattice[moving_particle+1] = 1
        # changing position Vector
        posR[moving_particle+1] = posR[moving_particle]+1
        posR[moving_particle] = 0
        # update state vector when particle jumps
        internal_state_vec[moving_particle+1] = internal_state_vec[moving_particle]
        internal_state_vec[moving_particle] = 0
        # change elongation_vector
        elongation_vector[moving_particle] = 0 #particle always leaves current position
        

    elseif moving_particle == 1 
        # compute the rightmost site we actually have
        front_end = min(l_ribosome + track_site, length(lattice))
        # only sum if that front_end ≥ 2 (because moving_particle+1 = 2)
        if 2 ≤ front_end && sum(lattice[2 : front_end]) == 0
            # now it’s safe to call the init function
            change_elong_vector_init_obc!(
                elongation_vector, lattice, moving_particle, rates, 
                l_ribosome, track_site, jammed, mobile, tot_part, num_clust
            )
        
            #hopping from site 1
            #change_elong_vector_init_obc(
            #    elongation_vector, lattice, moving_particle, rates, l_ribosome, track_site, 
            #    jammed, mobile, tot_part, num_clust
            #    )
            lattice[moving_particle+track_site] = 1
            # changing position Vector
            posR[moving_particle+track_site] = track_site + 1
            # update internal_state_vec
            internal_state_vec[moving_particle+track_site] = k₊
            internal_state_vec[moving_particle] = 0
            
            # function to change elongation_vector
            elongation_vector[moving_particle] = 0 #particle always leaves current position
        end

    elseif moving_particle == length(lattice)#- l_ribosome 
        # change lattice
        lattice[moving_particle] = 0
        # change position
        posR[moving_particle] = 0
        #end of cluster terminates
        if moving_particle+1 == c_end_pos_tracker[1] # plus 1 comes from the fact that we add +1 in the beginning of the function to the clsuter end tracker
            c_end_pos_tracker[1] = 0
            
        end
        
        if paused[1] == 0 && c_end_pos_tracker[1] == 0
            J_w[1] += 1
        elseif paused[1] >= 1 || c_end_pos_tracker[1] != 0
            J_c[1] += 1
        end

        J[1] += 1
        # update internal state vector
        mobile[1] -= 1
        internal_state_vec[moving_particle] = 0
        tot_part[1] -= 1
        
        # update elongation vector
        elongation_vector[moving_particle] = 0 #particle always leaves current position
        change_elong_vector_term_obc!(
            elongation_vector, moving_particle, internal_state_vec, 
            k₊, lattice, rates, l_ribosome, jammed, mobile, num_clust
            )
        
    end

    return moving_particle
end

function check_jammed(posR, lattice, l_ribosome, internal_state_vec, k₊)
    N = length(lattice)
    test = 0
    for site in posR
        if site == 0
            continue # skip the first site, it is always empty
        end
        ahead = site + l_ribosome
        if ahead ≤ N &&
           lattice[site] == 1 &&
           lattice[ahead] == 1 &&
           internal_state_vec[site] == k₊
            test += 1
        end
    end
    return test
end



function check_mobile(posR, lattice, l_ribosome, internal_state_vec, k₊)
    N = length(lattice)
    test = 0
    for site in posR
        ahead = site + l_ribosome
        if site == 0
            continue # skip the first site, it is always empty
        end
        if lattice[site] == 1 &&
           internal_state_vec[site] == k₊ &&
           (ahead > N || lattice[ahead] == 0)
            test += 1
        end
    end
    return test
end

function check_fct(check_jammed, lattice, posR, internal_state_vec, mobile, jammed, paused, l_ribosome, k₊)
    #check mobile + paused = total density
    #jammed = own fraction --> should always be less then mobile
    #jammed particles should be counted
    #elongation vector should be only occupied if mobile, non jammed, next site empty
    
    mistake = false
    if sum(lattice) != mobile[1] + paused[1] + jammed[1]
        mistake = true
    end
    if check_jammed(posR, lattice,l_ribosome, internal_state_vec, k₊) != jammed[1]
        mistake = true
    end
    if check_mobile(posR, lattice, l_ribosome, k₊) != mobile[1]
        mistake = true
    end
    return mistake
end

function correlation(lattice, l_ribosome)
    next_neighbor_occupancy = 0.0
    for site in eachindex(lattice[2:end-l_ribosome])
        if lattice[site] == 1 && lattice[site+l_ribosome] == 1
            next_neighbor_occupancy += 1
        end
    end
    avg_nn_correlation = next_neighbor_occupancy / length(lattice[2:end-l_ribosome])
    return avg_nn_correlation
end

function change_elong_vector_init_obc(
    elongation_vector, lattice, moving_particle, rates,
    l_ribosome, track_site, jammed, mobile, tot_part, num_clust)
    # no particle after jump
    if sum(lattice[1:1+l_ribosome+track_site]) == 0
        elongation_vector[moving_particle+track_site] = rates[moving_particle+track_site]
        
    else
        elongation_vector[moving_particle+track_site] = 0 
    end

    if lattice[1+track_site+l_ribosome] == 1
        jammed[1] += 1
    else
        mobile[1] += 1
        num_clust[1] += 1
    end

    # add 1 to particle number
    tot_part[1] += 1

end

function change_elong_vector_init_obc!(
    elongation_vector,
    lattice,
    moving_particle,   # always 1 when you call it
    rates,
    l_ribosome,
    track_site,
    jammed,
    mobile,
    tot_part,
    num_clust)
    N = length(lattice)

    # 1) Clamp the original slice 1 : (1 + l_ribosome + track_site) at N
    front_end = min(1 + l_ribosome + track_site, N)
    new_index = moving_particle + track_site   # = 1 + track_site

    if sum(lattice[1 : front_end]) == 0
        # footprint is clear
        if 1 ≤ new_index ≤ N
            elongation_vector[new_index] = rates[new_index]
        end
    else
        # footprint is blocked
        if 1 ≤ new_index ≤ N
            elongation_vector[new_index] = 0
        end
    end

    # 2) Decide jammed vs. mobile by checking exactly lattice[1 + track_site + l_ribosome], clamped
    check_index = 1 + track_site + l_ribosome
    if check_index ≤ N && lattice[check_index] == 1
        jammed[1] += 1
    else
        mobile[1]  += 1
        num_clust[1] += 1
    end

    # 3) Finally increment total particles
    tot_part[1] += 1
end


function change_elong_vector_bulk_obc(
    elongation_vector, moving_particle, lattice, 
    rates, internal_state_vec, k₊, l_ribosome, track_site, jammed, mobile, num_clust
    )

    if l_ribosome + track_site < moving_particle <= length(lattice)-1 #particle outside initiation region --> start looking "back"
        # particle moved
        if sum(lattice[moving_particle+1:moving_particle+l_ribosome+1])== 0
            elongation_vector[moving_particle+1] = rates[moving_particle+1] ## no particle in front
        end
        # check back if there is particle it becomes mobile
        if lattice[moving_particle-l_ribosome] == 1 && internal_state_vec[moving_particle-l_ribosome] == k₊ # particle behind must also be mobile
            elongation_vector[moving_particle-l_ribosome] = rates[moving_particle-l_ribosome]
            jammed[1] -= 1
            mobile[1] += 1
            
        end
        # if particle behind -> jump -> two clusters
        if lattice[moving_particle-l_ribosome] == 1
            num_clust[1] += 1
        end
        # check front
        if lattice[moving_particle+1+l_ribosome] == 0
            elongation_vector[moving_particle+1] = rates[moving_particle+1]
        else
            jammed[1] += 1
            mobile[1] -= 1
            num_clust[1] -= 1 #if particle infront after jump -> one cluster less
        end
    end
    # in initiation region
    if moving_particle <= l_ribosome + track_site
        
        if sum(lattice[moving_particle+1:moving_particle+l_ribosome+1]) == 0
            elongation_vector[moving_particle+1] = rates[moving_particle+1]
        end

        if moving_particle == l_ribosome + track_site # this means initation region is cleared
            elongation_vector[1] = rates[1]
        end

        if lattice[moving_particle+1+l_ribosome] == 1
            jammed[1] += 1
            mobile[1] -= 1
            num_clust[1] -= 1
        end
    end
end
function change_elong_vector_bulk_obc!(
    elongation_vector,
    moving_particle,
    lattice,
    rates,
    internal_state_vec,
    k₊,
    l_ribosome,
    track_site,
    jammed,
    mobile,
    num_clust)
    N = length(lattice)

    # === CASE A: particle is outside the initiation region ===
    if l_ribosome + track_site < moving_particle ≤ N - 1
        #
        # 1) “Look ahead” exactly l_ribosome+1 sites, but clamp end to N:
        front_start = moving_particle + 1
        front_end   = min(moving_particle + l_ribosome + 1, N)

        if front_start ≤ N
            # If there are no particles in [front_start:front_end],
            # we can set elongation_vector at (moving_particle+1).
            if sum(lattice[front_start : front_end]) == 0
                if moving_particle + 1 ≤ N
                    elongation_vector[moving_particle + 1] = rates[moving_particle + 1]
                end
            end
        end

        #
        # 2) “Look behind” exactly at (moving_particle - l_ribosome), if ≥ 1:
        back_index = moving_particle - l_ribosome
        if back_index ≥ 1
            # If a ribosome is exactly l_ribosome behind AND is in state k₊, it becomes mobile.
            if lattice[back_index] == 1 && internal_state_vec[back_index] == k₊
                elongation_vector[back_index] = rates[back_index]
                jammed[1]  -= 1
                mobile[1]  += 1
            end

            # If there is any ribosome at back_index → that created two clusters
            if lattice[back_index] == 1
                num_clust[1] += 1
            end
        end

        #
        # 3) Re‐check “just in front of the new position” (i.e. index = moving_particle+1+l_ribosome)
        front_block_index = moving_particle + 1 + l_ribosome
        if front_block_index ≤ N
            if lattice[front_block_index] == 0
                # No ribosome blocking—ensure moved particle (at moving_particle+1) is mobile:
                elongation_vector[moving_particle + 1] = rates[moving_particle + 1]
            else
                # There is a ribosome exactly l_ribosome ahead → moved particle is jammed
                jammed[1]    += 1
                mobile[1]    -= 1
                num_clust[1] -= 1
            end
        end
    end


    # === CASE B: particle is still within the initiation region ===
    if moving_particle ≤ l_ribosome + track_site
        #
        # 1) “Look ahead” exactly l_ribosome+1 sites, clamp to N:
        front_start = moving_particle + 1
        front_end   = min(moving_particle + l_ribosome + 1, N)

        if front_start ≤ N && sum(lattice[front_start : front_end]) == 0
            elongation_vector[moving_particle + 1] = rates[moving_particle + 1]
        end

        #
        # 2) If moving_particle was exactly at the last site of the initiation region,
        #    that region has now cleared, so initiation can fire at site 1:
        if moving_particle == l_ribosome + track_site
            elongation_vector[1] = rates[1]
        end

        #
        # 3) Check “just in front” at index = moving_particle + 1 + l_ribosome
        check_index = moving_particle + 1 + l_ribosome
        if check_index ≤ N && lattice[check_index] == 1
            jammed[1]    += 1
            mobile[1]    -= 1
            num_clust[1] -= 1
        end
    end
end

function change_elong_vector_term_obc!(
    elongation_vector,
    moving_particle,
    internal_state_vec,
    k₊,
    lattice,
    rates,
    l_ribosome,
    jammed,
    mobile,
    num_clust)
    N = length(lattice)
    back_index = moving_particle - l_ribosome

    # Only look “behind” if back_index ≥ 1
    if back_index ≥ 1
        if lattice[back_index] == 1 && internal_state_vec[back_index] == k₊
            # the ribosome behind is mobile → it becomes free to move
            elongation_vector[back_index] = rates[back_index]
            jammed[1]  -= 1
            mobile[1]  += 1

        elseif lattice[back_index] == 0
            # nothing directly behind → one cluster ends
            num_clust[1] -= 1
        end
    end
end

function change_elong_vector_term_obc(
    elongation_vector, moving_particle, 
    internal_state_vec, k₊, lattice, rates,
    l_ribosome, jammed, mobile, num_clust
    )
    # if particle behind terminated particle is not paused it becomes mobile
    if lattice[moving_particle-l_ribosome] == 1 && internal_state_vec[moving_particle-l_ribosome] == k₊
        elongation_vector[moving_particle-l_ribosome] = rates[moving_particle-l_ribosome]
        jammed[1] -= 1
        mobile[1] += 1
    
    elseif lattice[moving_particle-l_ribosome] == 0
        num_clust[1] -= 1

    end

end

function change_internal_state_obc!(
    internal_state_vec,
    w_state,
    posR,
    lattice,
    rates,
    elongation_vector,
    k₋,
    k₊,
    mobile,
    paused,
    jammed,
    l_ribosome,
    c_end_pos_tracker,
    pos_first_paused_particle,
    delta_t)

    N = length(lattice)
    # Pick a particle to switch:
    switching_pos = sample(posR, w_state)

    # Compute the index that is "l_ribosome" ahead of switching_pos:
    ahead_index = switching_pos + l_ribosome

    if internal_state_vec[switching_pos] == k₋
        # ------------------ PAUSED → MOBILE/JAMMED ------------------
        internal_state_vec[switching_pos] = k₊
        paused[1] -= 1

        # If we just went from 1 paused → 0 paused, find the cluster end:
        if paused[1] == 0
            # find_cluster_end should set c_end_pos_tracker[1] based on switching_pos
            find_cluster_end(switching_pos, posR, l_ribosome, c_end_pos_tracker)
        end

        # Decide if the “now‐unpaused” particle becomes jammed or mobile:
        if ahead_index ≤ N && lattice[ahead_index] == 1
            # There's a ribosome exactly l_ribosome ahead, so the unpaused particle is jammed
            jammed[1] += 1
        else
            # Either (ahead_index > N) = off the lattice, or it is empty on the lattice
            mobile[1] += 1
        end

        # If this particle is now mobile _and_ it’s not past the last site, it contributes to elongation:
        if (switching_pos ≤ N) && (ahead_index > N || lattice[ahead_index] == 0)
            # when ahead_index > N, we treat “off‐end” as empty—so it’s free to go
            elongation_vector[switching_pos] = rates[switching_pos]
        end

    else
        # --------------- MOBILE/JAMMED → PAUSED ------------------
        internal_state_vec[switching_pos] = k₋
        elongation_vector[switching_pos] = 0
        paused[1] += 1

        if paused[1] == 1
            # First paused in this cluster: record time
            pos_first_paused_particle[1] += delta_t
            c_end_pos_tracker[1] = 0
        end

        # If there was a ribosome l_ribosome ahead, the switching particle loses its jammed status if not it loses its mobile status:
        if ahead_index ≤ N && lattice[ahead_index] == 1
            #was jammed becomes paused
            jammed[1] -= 1
        else
            # Either “off‐end” (ahead_index > N) or no neighbor → so we removed one mobile
            mobile[1] -= 1
        end
    end

    return switching_pos
end


function change_internal_state_obc(
    internal_state_vec, posR, lattice, rates, elongation_vector, 
    k₋, k₊, mobile, paused, jammed, l_ribosome,c_end_pos_tracker, pos_first_paused_particle, delta_t)

    switching_pos = sample(posR, Weights(internal_state_vec))

    if internal_state_vec[switching_pos] == k₋ # switch from paused to mobile or jammed
        
        internal_state_vec[switching_pos] = k₊ # the particle at said position is non paused now and can become paused with a rate f

        paused[1] -= 1 #substract paused particle due to switching

        #if paused[1] == 1 # means we are down to the last piece of the cluster
        #    # need a function here that counts till the end of the cluster
        #    c_end_pos_tracker[1] = switching_pos # the end of the cluster 
        #end

        if paused[1] == 0 # means we are down to the last piece of the cluster
            # need a function here that counts till the end of the cluster
            #one_to_zero[1] += 1
            find_cluster_end(switching_pos, posR,l_ribosome,c_end_pos_tracker)
        end

        if lattice[switching_pos+l_ribosome] == 1 # if neighbor site of previous paused particle is occupied a jammed particle is added
            jammed[1] += 1  # add a jammed particle
        else   # if the neighboring site is unoccupied the particle has to be mobile 
            mobile[1] += 1 # add mobile particle
        end

        if switching_pos != length(lattice) && lattice[switching_pos+l_ribosome] == 0  # excluding check for sites beyond lattice (due to addtion of larger particles)
            elongation_vector[switching_pos] = rates[switching_pos] # particle in bulk contributes
        elseif switching_pos == length(lattice) #mean last site is choosen & occupied --> particle from paused to mobile
            elongation_vector[switching_pos] = rates[switching_pos] # the particle becomes active and contributes to elongation
        end

    else # means a mobile or jammed particle state was chosen which has to switch to a paused state
    
        internal_state_vec[switching_pos] = k₋ # particle is paused now with a life time of τ
        elongation_vector[switching_pos] = 0 # paused particle does not contribute to elongation vector

        paused[1] += 1 # adding to paused state
        
        if paused[1] == 1 # we switch from zero to one paused particle
            # need a function here that counts till the end of the cluster
            #if c_end_pos_tracker[1] != 0.0
            #    zero_to_one[1] += 1
            #end
            #pos_first_paused_particle[1] += delta_t 
            c_end_pos_tracker[1] = 0 # this means we added a paused particle to the system we are no longer in the clogged state
            
        end

        if lattice[switching_pos+l_ribosome] == 1 #next site is occuppied
            jammed[1] -= 1 # add jammed
        else # next site is empty
            mobile[1] -= 1 # add mobile
        end
    end

    if paused[1] == 1
        c_end_pos_tracker[1] = switching_pos
    end

    return switching_pos
end

function find_cluster_end(switching_pos, posR,l_ribosome,c_end_pos_tracker)
    
    end_not_found = true
    i = 0
    while end_not_found
        if switching_pos-l_ribosome*i >= 2 && posR[switching_pos-l_ribosome*i] != 0.0 
            i += 1
        else
            end_not_found = false
        end
    end
    c_end_pos_tracker[1] = switching_pos - l_ribosome*(i-1)
end

function reset_observables(lattice, rates, elongation_vector, internal_state_vec, posR, tot_part, mobile, jammed, paused)
    lattice .= 0
    elongation_vector .= 0
    elongation_vector[1] = rates[1]
    internal_state_vec .= 0
    posR .= 0
    posR[1] = 1
    tot_part .= 0
    mobile .= 0
    jammed .= 0
    paused .= 0
end

function Gillespie_obc(L, l_ribosome, track_site, deg_t, deg_g, init, term, elong, k₋, k₊, run_time, starting_t, delta_t; kymo=false)
    
    if k₋ == k₊
        k₊ *= 0.99999999
    end

    number_of_measurements = Int64((run_time-starting_t) / delta_t) # number of measurements

    # number of different particle classes
    mobile = [0.0] 
    paused = [0.0]
    jammed = [0.0]
    
    tot_part = [0.0]
    tot_weighted = [0.0]
    tot_weighted_p = [0.0]
    tot_weighted_unp = [0.0]

    tot_vec_p = Vector{Float64}(undef,number_of_measurements)
    tot_vec_unp = Vector{Float64}(undef,number_of_measurements)
    tot_vec = Vector{Float64}(undef,number_of_measurements)

    num_clust = [0.0]

    unpausing_counter = [0.0]
    unpausing_rate  = [0.0]

    pausing_counter = [0.0]
    pausing_rate_1 = [0.0]
    
    t_c_vec = Vector{Float64}(undef, number_of_measurements)
    t_w_vec = Vector{Float64}(undef, number_of_measurements)  
    time_vec = Vector{Float64}(undef, number_of_measurements)

    # ρ_vec = Vector{Float64}(undef,number_of_measurements) # density vector
    # ρ_weighted = [0.0] # density after weighted by time step the configuration

    mobile_vec = Vector{Float64}(undef, number_of_measurements) # mobile vector
    mobile_weighted = [0.0] # mobile fraction of particles

    jammed_vec= Vector{Float64}(undef,number_of_measurements)# jammed vector
    jammed_weighted = [0.0]# jammed fraction of particles

    paused_vec = Vector{Float64}(undef,number_of_measurements)# paused vector
    paused_weighted = [0.0]# paused fraction of particles

    current = Vector{Float64}(undef,number_of_measurements) # current vector
    w_current = Vector{Float64}(undef,number_of_measurements)
    c_current = Vector{Float64}(undef,number_of_measurements)

    J_w = [0] # TASEP current
    J_c = [0] # pTASEP current
    J = [0] # current counter at last site
    J_c_end = [0]
    c_end_pos_tracker = [0.0]



    lattice = create_lattice_l(L,l_ribosome) # create lattice (lattice= start site(always empty) + init_region (length(l_ribosome)) + L + l_ribosome)
    rates = create_rates_l(init, term, elong, L) # rate vector
    posR = positionRibosomes(lattice) # position vector
    elongation_vector = get_elongation_vector_obc(lattice, rates) # which particles contribute to configuration
    internal_state_vec = get_internal_state_vec_obc(lattice, k₊, mobile) # which particle contributes to configuration via its internal state

    w_elong  = Weights(elongation_vector)
    w_elong.values === elongation_vector
    
    w_state  = Weights(internal_state_vec)
    w_state.values === internal_state_vec

    cluster_num_bucket_p = initiate_cluster_buckets_obc(lattice)
    cluster_num_bucket_unp = initiate_cluster_buckets_obc(lattice)
    cluster_len_bucket = initiate_cluster_buckets_obc(lattice)
    paused_dist_bucket = initiate_cluster_buckets_obc(lattice)
    pos_first_paused_particle_bucket = initiate_cluster_buckets_obc(lattice)
    #time = [] # track how many steps in one delta_t
    #time_temp = [] 
    t = 0.0
    #ij = 0
    while t < starting_t # if time is below starting_t, no measurements are taken
        #ij += 1
        #w_elong.values .= elongation_vector # sum of all the rates in the elongation vector
        w_elong.sum    = sum(elongation_vector)

        #w_state.values .= internal_state_vec
        w_state.sum    = sum(internal_state_vec)

    
        
        deg_rate = deg_g
        #if paused[1] == 0 #no paused particle the degrdation is equal to the general one
        #    deg_rate = deg_g 
        #else # at least a single paused particle the degradation rate is targeted
        #    deg_rate = deg_t
        #end

        total_sum = w_elong.sum  + w_state.sum  + deg_rate
        RV = rand() # Random variable to choose between elongation and internal state switching

        if RV <= w_elong.sum /total_sum # elongation is chosen
            elongation_process_obc(
                elongation_vector, w_elong, internal_state_vec, posR, lattice, rates, 
                J, k₊, mobile, l_ribosome, track_site, jammed, tot_part,num_clust,
                J_c, J_w, paused, c_end_pos_tracker
                )

        elseif w_elong.sum /total_sum < RV <= (w_elong.sum  + w_state.sum) / total_sum  # switching happens
            change_internal_state_obc!(
                internal_state_vec, w_state, posR, lattice, rates, elongation_vector, 
                k₋, k₊, mobile, paused, jammed, l_ribosome,c_end_pos_tracker,pos_first_paused_particle_bucket,delta_t)
        else
            reset_observables(lattice, rates, elongation_vector, internal_state_vec, posR, tot_part, mobile, jammed, paused)
                
        end
        time_added = rand(Exponential(1/total_sum))
        t += time_added # add time depending on the lattice configuration
        #print("$t\n")
        #print("$ij\n")
    end

    J[1] = 0 # reset current
    J_w[1] = 0
    J_c[1] = 0

    deg_counter_paused = [0.0]
    deg_counter_non_paused = [0.0]

    pausing_counter[1] = 0.0 # reset pausing counter
    pausing_rate_1[1] = 0.0 # reset rate

    unpausing_counter[1] = 0.0 # reset unpausing counter
    unpausing_rate[1] = 0.0 # reset rate

    c_end_pos_tracker[1] = 0.0
    J_c_end[1] = 0.0

    t_c = 0.0
    t_w = 0.0
    t_c_old = 0.0
    t_w_old = 0.0
    
    t = 0.0 # reset time counter

    for i in 1:number_of_measurements
        #j = 1
        #print("$j")
        while t <= delta_t # if time is below delta_t (time interval)

            #w_elong.values .= elongation_vector # sum of all the rates in the elongation vector
            w_elong.sum    = sum(elongation_vector)

            #w_state.values .= internal_state_vec
            w_state.sum    = sum(internal_state_vec)


            deg_rate = deg_g



            total_sum = w_elong.sum  + w_state.sum + deg_rate
            RV = rand() # Random variable to choose between elongation and internal state switching
            time_added = rand(Exponential(1/total_sum)) # observeables are weigthed by their configuration
            
            if paused[1] == 0 && c_end_pos_tracker[1] == 0

                t_w += time_added
                tot_weighted_unp[1] += tot_part[1]/length(posR) * time_added

            elseif paused[1] >= 1 || c_end_pos_tracker[1] != 0

                t_c += time_added
                tot_weighted_p[1] += tot_part[1]/length(posR) * time_added

            end #no paused particle add time to clogged 
                

            t += time_added # add to time
            
            mobile_weighted[1] += mobile[1]/length(posR) * time_added # fraction of mobile particles * time from previous configuration
            paused_weighted[1] += paused[1]/length(posR) * time_added # fractopm of paused particles
            jammed_weighted[1] += jammed[1]/length(posR) * time_added # fraction of jammed particles
            tot_weighted[1] += tot_part[1]/length(posR) * time_added

            #if sum(lattice) != tot_part[1]
            #    print("$(sum(lattice)) = $(tot_part[1])\n")
            #    print("holy shit I am fucked\n")
            #end
            #cluster_num_bucket[num_clust[1]] += time_added
            #check_m = check_mobile(posR, lattice, l_ribosome, internal_state_vec, k₊)
            #check_j = check_jammed(posR, lattice, l_ribosome, internal_state_vec, k₊)
            #if check_j != jammed[1] || check_m != mobile[1]
            #    print("Houston we have a problem")
            #    break
            #end

            #find_all_cluster_lengths_obc(
            #    lattice, 
            #    cluster_len_bucket, 
            #    cluster_num_bucket_unp,
            #    cluster_num_bucket_p, 
            #    time_added, 
            #    num_clust,
            #    paused,
            #    l_ribosome,
            #    c_end_pos_tracker
            #    )
            
            #paused_dist_bucket[paused[1]] += time_added
            
            #cluster_len_num(time_added, cluster_num_bucket, cluster_len_bucket, lattice, num_clust)
            #sort_clusters_length_and_dist(lattice, paused, internal_state_vec, paused_dist_bucket, k₊, time_added, k₋)
            # ρ_weighted[1] += mean(lattice[2:end-l_ribosome]) * time_added # particle density

           
            if RV <= w_elong.sum /total_sum # elongation is chosen
                elongation_process_obc(
                    elongation_vector, w_elong, internal_state_vec, posR, lattice, rates, 
                    J, k₊, mobile, l_ribosome, track_site, jammed, tot_part,num_clust,
                    J_c, J_w, paused,c_end_pos_tracker
                    )

            elseif w_elong.sum /total_sum < RV <= (w_elong.sum  + w_state.sum) / total_sum  # switching happens
                change_internal_state_obc!(
                    internal_state_vec, w_state, posR, lattice, rates, elongation_vector, 
                    k₋, k₊, mobile, paused, jammed, l_ribosome,c_end_pos_tracker,pos_first_paused_particle_bucket,delta_t)
            else

                if paused[1] == 0 
                    deg_counter_non_paused[1] += 1
                else
                    deg_counter_paused[1] += 1
                end

                reset_observables(lattice, rates, elongation_vector, internal_state_vec, posR, tot_part, mobile, jammed, paused)
                
            end
            
            
            # nn_occupancy[1] += correlation(lattice, l_ribosome) * time_added
            #print("J: $(J[1])\n")
            #j += 1
        end
        
        # divide the measurements by delta_t/batch size
        current[i]= J[1]/t 
        w_current[i] = J_w[1] / t
        c_current[i] = J_c[1] / t

        # ρ_vec[i] = ρ_weighted[1]/delta_t
        
        mobile_vec[i] = mobile_weighted[1]/t
        paused_vec[i] = paused_weighted[1]/t
        jammed_vec[i] = jammed_weighted[1]/t
        tot_vec[i] = tot_weighted[1]/t
        tot_vec_unp[i] = tot_weighted_unp[1] /t
        tot_vec_p[i] = tot_weighted_p[1] / t
        # correlation_vec[i] = nn_occupancy[1]/delta_t - (ρ_vec[i])^2

        time_vec[i] = t
        t = 0
        t_w_vec[i] = t_w 
        t_w = 0
        t_c_vec[i] = t_c 
        t_c = 0


        #reset the counters so they are ready for new delta_t
        J[1] = 0.0
        J_w[1] = 0.0
        J_c[1] = 0.0

        # ρ_weighted[1] = 0.0
        mobile_weighted[1] = 0.0
        paused_weighted[1] = 0.0
        jammed_weighted[1] = 0.0

        tot_weighted[1] = 0.0
        tot_weighted_p[1] = 0.0
        tot_weighted_unp[1] = 0.0
        # nn_occupancy[1] = 0.0
    
        
    end

    # take the mean off all the time batches
    J_mean::Float64 = mean(current)
    J_w_mean::Float64 = mean(w_current)
    J_c_mean::Float64 = mean(c_current)

    ρ_vec::Vector{Float64} = vcat(
        mean(mobile_vec), 
        mean(paused_vec), 
        mean(jammed_vec), 
        mean(tot_vec),  
        mean(tot_vec_unp),
        mean(tot_vec_p)
    ) 

    #P_res = unpausing_counter[1] / pausing_counter[1]
    
    P_paused_deg = deg_counter_paused[1] / (deg_counter_non_paused[1] + deg_counter_paused[1])

    #pausing_r = pausing_rate_1[1] #/ delta_t
    #unpausing_r = unpausing_rate[1]# /delta_t

    #time = hcat(mean(time_vec), mean(t_w_vec), mean(t_c_vec))
    time_in_w = mean(t_w_vec)
    time_in_p = mean(t_c_vec)
    total_time = mean(time_vec)

    paused_particle_distribution = Dict(key=>value/total_time for (key, value) in paused_dist_bucket)
    cluster_length_distribution = Dict(key=>value/total_time for (key, value) in cluster_len_bucket)
    cluster_number_distribution_unp = Dict(key=>value/total_time for (key, value) in cluster_num_bucket_unp)
    cluster_number_distribution_p = Dict(key=>value/total_time for (key, value) in cluster_num_bucket_p)

    # corr_mean::Float64 = mean(correlation_vec)
    #clean_clusters(cluster_num_bucket)
    #clean_clusters(cluster_len_bucket)
    #clean_clusters(paused_dist_bucket)
    #clusters = [cluster_num_bucket cluster_len_bucket paused_dist_bucket]
     
    #if kymo == true
    #    return J_mean, ρ_matrix, t_vec, kymo_matr, intern_matr
    #end
    
    #return J_mean, ρ_matrix, total_time
    return [J_w_mean, J_c_mean, J_mean], ρ_vec, paused_particle_distribution, cluster_length_distribution , cluster_number_distribution_unp, cluster_number_distribution_p, [time_in_w, time_in_p, total_time], P_paused_deg 
end

function Gillespie_obc_kymo_steps(L, l_ribosome, init, term, elong, k₋, k₊, total_steps, starting_steps, record_interval)
    if k₋ == k₊
        k₊ *= 0.99999999
    end

    mobile = [0.0] 
    paused = [0.0]
    jammed = [0.0]

    # Initialize the lattice and rates
    lattice = create_lattice_l(L,l_ribosome) # create lattice (lattice= start site(always empty) + init_region (length(l_ribosome)) + L + l_ribosome)
    rates = create_rates_l(init, term, elong, L) # rate vector
    posR = positionRibosomes(lattice) # position vector
    elongation_vector = get_elongation_vector_obc(lattice, rates) # which particles contribute to configuration
    internal_state_vec = get_internal_state_vec_obc(lattice, k₊, mobile) # which particle contributes to configuration via its internal state

    w_elong  = Weights(elongation_vector)
    w_elong.values === elongation_vector
    
    w_state  = Weights(internal_state_vec)
    w_state.values === internal_state_vec

    # Initialize step counters and kymograph storage
    step = 0
    kymo_matr = lattice
    intern_matr = internal_state_vec
    steps_vec = []

    # Burn-in period (no measurements taken)
    while step < starting_steps
        elong_sum = sum(elongation_vector)
        state_sum = sum(internal_state_vec)

        total_sum = elong_sum + state_sum 
        RV = rand()
        if RV <= elong_sum / total_sum
            elongation_process_obc(
                    elongation_vector, w_elong, internal_state_vec, posR, lattice, rates, 
                    [0.0], k₊, [0.0], l_ribosome, track_site, [0.0], [0.0],[0.0],
                    [0.0], [0.0], [0.0],[0]
                    )
        elseif RV <= (elong_sum + state_sum) / total_sum
             change_internal_state_obc!(
                internal_state_vec, w_state, posR, lattice, rates, elongation_vector, 
                k₋, k₊, [0.0], [0.0], [0.0], l_ribosome,[0],[0], 1)

        end
        step += 1
    end

    # Simulation with kymograph recording
    step = 0
    total_time = 0.0  # Optional: to keep track of time
    time_vec = []

    while step < total_steps
        elong_sum = sum(elongation_vector)
        state_sum = sum(internal_state_vec)

        total_sum = elong_sum + state_sum
        RV = rand()
        time_added = rand(Exponential(1 / total_sum))
        total_time += time_added  # Optional: accumulate time

        # Record data at specified steps
        if step % record_interval == 0
            push!(steps_vec, step)
            push!(time_vec, total_time)  # Optional: store corresponding times
            kymo_matr = hcat(kymo_matr, lattice)
            intern_matr = hcat(intern_matr, internal_state_vec)
        end

        if RV <= elong_sum / total_sum
            elongation_process_obc(
                elongation_vector, w_elong, internal_state_vec, posR, lattice, rates, 
                [0.0], k₊, [0.0], l_ribosome, track_site, [0.0], [0.0],[0.0],
                [0.0], [0.0], [0.0],[0])
        elseif RV <= (elong_sum + state_sum) / total_sum
             change_internal_state_obc!(
                internal_state_vec, w_state, posR, lattice, rates, elongation_vector, 
                k₋, k₊, [0.0], [0.0], [0.0], l_ribosome,[0],[0], time_added)
        end
        step += 1
    end

    return steps_vec, time_vec, kymo_matr, intern_matr
end

function Gillespie_obc_kymo_time(L, l_ribosome, init, term, elong, k₋, k₊, run_time, starting_t, delta_t)
    if k₋ == k₊
        k₊ *= 0.99999999
    end

    mobile = [0.0] 
    paused = [0.0]
    jammed = [0.0]

    # Initialize the lattice and rates
    lattice = create_lattice_l(L,l_ribosome) # create lattice (lattice= start site(always empty) + init_region (length(l_ribosome)) + L + l_ribosome)
    rates = create_rates_l(init, term, elong, L) # rate vector
    posR = positionRibosomes(lattice) # position vector
    elongation_vector = get_elongation_vector_obc(lattice, rates) # which particles contribute to configuration
    internal_state_vec = get_internal_state_vec_obc(lattice, k₊, mobile) # which particle contributes to configuration via its internal state

    w_elong  = Weights(elongation_vector)
    w_elong.values === elongation_vector
    
    w_state  = Weights(internal_state_vec)
    w_state.values === internal_state_vec

    # Initialize time and kymograph storage
    t = 0.0
    kymo_matr = lattice
    intern_matr = internal_state_vec
    t_vec = []

    # Burn-in period (no measurements taken)
    while t < starting_t
        w_elong.sum    = sum(elongation_vector)

        #w_state.values .= internal_state_vec
        w_state.sum    = sum(internal_state_vec)


        total_sum = w_elong.sum  + w_state.sum 
        RV = rand()
        if RV <= w_elong.sum / total_sum
            elongation_process_obc(
                    elongation_vector, w_elong, internal_state_vec, posR, lattice, rates, 
                    [0.0], k₊, [0.0], l_ribosome, track_site, [0.0], [0.0],[0.0],
                    [0.0], [0.0], [0.0],[0]
                    )

        elseif RV <= (w_elong.sum + w_state.sum ) / total_sum
             change_internal_state_obc!(
                internal_state_vec, w_state, posR, lattice, rates, elongation_vector, 
                k₋, k₊, [0.0], [0.0], [0.0], l_ribosome,[0.0],[0.0], 0.0)
        end
        time_added = rand(Exponential(1 / total_sum))
        t += time_added
    end

    # Simulation with kymograph recording
    t = 0.0
    t_step = delta_t  # Adjust the time step as needed
    next_record_time = t_step

    while t < run_time
        w_elong.sum    = sum(elongation_vector)

        #w_state.values .= internal_state_vec
        w_state.sum    = sum(internal_state_vec)


        total_sum = w_elong.sum + w_state.sum 
        RV = rand()
        time_added = rand(Exponential(1 / total_sum))
        t += time_added

        if t >= next_record_time
            push!(t_vec, t)
            kymo_matr = hcat(kymo_matr, lattice)
            intern_matr = hcat(intern_matr, internal_state_vec)
            next_record_time += t_step
        end

        if RV <= w_elong.sum  / total_sum
                elongation_process_obc(
                    elongation_vector, w_elong, internal_state_vec, posR, lattice, rates, 
                    [0.0], k₊, [0.0], l_ribosome, track_site, [0.0], [0.0],[0.0],
                    [0.0], [0.0], [0.0],[0]
                    )
        elseif RV <= (w_elong.sum   + w_state.sum ) / total_sum
             change_internal_state_obc!(
                internal_state_vec, w_state, posR, lattice, rates, elongation_vector, 
                k₋, k₊, [0.0], [0.0], [0.0], l_ribosome,[0.0],[0.0], time_added)
        end
    end

    return t_vec, kymo_matr, intern_matr
end

function find_all_cluster_lengths_obc(
    vector, 
    cluster_lengths, 
    cluster_numbers_unp, 
    cluster_numbers_p,
    time_added, 
    num_clust,
    paused,
    ribosome_size,  # Number of sites occupied by a ribosome (e.g., 10 codons)
    c_end_pos_tracker
    )

    current_length = 0
    cluster_number = 0

    if sum(vector) == 0
        cluster_lengths[0] += time_added
    else
        i = 1
        while i <= length(vector)
            if vector[i] == 1
                current_length += 1  # particle
                i += ribosome_size  # Skip `ribosome_size` positions ahead
            elseif current_length > 0
                # Register cluster length
                cluster_lengths[current_length] += time_added / num_clust[1]
                cluster_number += 1

                # Reset current length for next cluster
                current_length = 0
                i += 1
            else
                i += 1
            end
        end

        # If there was a cluster at the end of the vector
        if current_length > 0
            cluster_lengths[current_length] += time_added / num_clust[1]
            cluster_number += 1
        end

        # Update the number of clusters based on pause state

        

    end
    if paused[1] == 0 && c_end_pos_tracker[1] == 0
        cluster_numbers_unp[cluster_number] += time_added
    elseif paused[1] >= 1 || c_end_pos_tracker[1] != 0
        cluster_numbers_p[cluster_number] += time_added
    end
    # Optionally: sanity check to ensure cluster count matches expectations
    if cluster_number != num_clust[1]
        print("Warning: Something is wrong with the cluster count")
    end
end

function run_until_first_pause_steady_state(
    l_ribosome::Int, k_minus::Float64, k_plus::Float64,
    lattice::Vector{Int}, rates::Vector{Float64}, posR::Vector{Int}, elongation_vector::Vector{Float64},
    mobile::Vector{Float64}, jammed::Vector{Float64}, paused::Vector{Float64}, internal_state_vec::Vector{Float64}, track_site::Any,
    pos_first_paused_particle_bucket::Vector{Float64}, tot_part::Vector{Float64}, num_clust::Vector{Float64}, 
    J::Vector{Int}, J_c::Vector{Int}, J_w::Vector{Int}, c_end_pos_tracker::Vector{Float64})::Int
    
    @inbounds begin
        # Precompute the stopping condition once:
        max_time = 10 * length(posR) / rates[end]
        t = 0.0
        
        while t < max_time
            elong_sum = 0.0
            @inbounds @simd for i in eachindex(elongation_vector)
                elong_sum += elongation_vector[i]
            end
            state_sum = 0.0
            @inbounds @simd for i in eachindex(internal_state_vec)
                state_sum += internal_state_vec[i]
            end
            total_sum = elong_sum + state_sum

            # Draw time increment from exponential:
            # Instead of rand(Exponential(1/total_sum)), we can do:
            # time_added = randexp() / total_sum
            # where randexp() = -log(rand())
            time_added = -log(rand())/total_sum
            t += time_added

            RV = rand()
            if RV <= elong_sum/total_sum
                elongation_process_obc(
                    elongation_vector, internal_state_vec, posR, lattice, rates, 
                    J, k_plus, mobile, l_ribosome, track_site, jammed, tot_part,num_clust,
                    J_c, J_w, paused, c_end_pos_tracker
                )
            end
        end

        # After the time limit
        while true
            elong_sum = 0.0
            @inbounds @simd for i in eachindex(elongation_vector)
                elong_sum += elongation_vector[i]
            end
            state_sum = 0.0
            @inbounds @simd for i in eachindex(internal_state_vec)
                state_sum += internal_state_vec[i]
            end
            total_sum = elong_sum + state_sum

            time_added = -log(rand())/total_sum
            RV = rand()

            if RV <= elong_sum/total_sum
                elongation_process_obc(
                    elongation_vector, internal_state_vec, posR, lattice, rates, 
                    J, k_plus, mobile, l_ribosome, track_site, jammed, tot_part,num_clust,
                    J_c, J_w, paused, c_end_pos_tracker
                )
            else
                switching_pos = change_internal_state_obc(
                    internal_state_vec, posR, lattice, rates, elongation_vector, 
                    k_minus, k_plus, mobile, paused, jammed, l_ribosome, c_end_pos_tracker, 
                    pos_first_paused_particle_bucket, time_added
                )
                return switching_pos
            end
        end
    end
end

function density_profile(
    l_ribosome::Int, k_plus::Float64,
    lattice::Vector{Int}, rates::Vector{Float64}, 
    posR::Vector{Int}, elongation_vector::Vector{Float64},
    mobile::Vector{Float64}, jammed::Vector{Float64}, 
    paused::Vector{Float64}, internal_state_vec::Vector{Float64}, track_site::Any, 
    tot_part::Vector{Float64}, num_clust::Vector{Float64}, 
    J::Vector{Int}, J_c::Vector{Int}, J_w::Vector{Int}, 
    c_end_pos_tracker::Vector{Float64}, number_of_steps::Int64
    )
    
    t = 0
    dt_vec = Vector{Float64}(undef, number_of_steps)
    total_time = Vector{Float64}(undef,number_of_steps)
    density_profile = Matrix{Int}(undef, number_of_steps, length(lattice))

    # After the time limit
    #for step in 1:number_of_steps
    while J[1] == 0
        elong_sum = sum(elongation_vector)
        total_sum = elong_sum 
        dt = rand(Exponential(1/total_sum))

        #dt_vec[step] = dt
        if sum(lattice) != 0
            t += dt
        end
        #total_time[step] = t

        #density_profile[step,1:end] = lattice

        elongation_process_obc(
            elongation_vector, internal_state_vec, posR, lattice, rates, 
            J, k_plus, mobile, l_ribosome, track_site, jammed, tot_part,num_clust,
            J_c, J_w, paused, c_end_pos_tracker
        )
        

    end
    #final_data = hcat(dt_vec, total_time, Float64.(density_profile))
    #return final_data
    return t
end

function gather_density_profile(
    L::Int, l_ribosome::Int, track_site, init::Float64, term::Float64, elong::Float64, 
    k_minus::Float64, k_plus::Float64, number_of_steps::Int64)
    # Preallocate arrays once:
    lattice = create_lattice_l(L, l_ribosome)             # ensure these return arrays with concrete types
    rates = create_rates_l(init, term, elong, L, l_ribosome)
    posR = positionRibosomes(lattice, l_ribosome)
    elongation_vector = get_elongation_vector_obc(lattice, rates, l_ribosome)

    mobile = [0.0]
    paused = [0.0]
    jammed = [0.0]
    internal_state_vec = get_internal_state_vec_obc(lattice, k_plus, mobile, l_ribosome)

    pos_first_paused_particle_bucket = [0.0]
    tot_part = [0.0]
    num_clust = [0.0]
    J_c = [0]
    J_w = [0]
    J = [0]
    c_end_pos_tracker = [0.0]
    

    profile = density_profile(
        l_ribosome, k_plus, lattice, rates, posR, elongation_vector,
        mobile, jammed, paused, internal_state_vec, track_site, tot_part,
        num_clust, J, J_c, J_w, c_end_pos_tracker, number_of_steps
    )

    return profile
end

function run_until_first_pause_profile(
    l_ribosome::Int, k_minus::Float64, k_plus::Float64,
    lattice::Vector{Int}, rates::Vector{Float64}, posR::Vector{Int}, elongation_vector::Vector{Float64},
    mobile::Vector{Float64}, jammed::Vector{Float64}, paused::Vector{Float64}, internal_state_vec::Vector{Float64}, track_site::Any,
    pos_first_paused_particle_bucket::Vector{Float64}, tot_part::Vector{Float64}, num_clust::Vector{Float64}, 
    J::Vector{Int}, J_c::Vector{Int}, J_w::Vector{Int}, c_end_pos_tracker::Vector{Float64})
    
    t = 0
    density_profile = lattice
    @inbounds begin
        
        # After the time limit
        while true
        
            elong_sum = sum(elongation_vector)
            state_sum = sum(internal_state_vec)
            total_sum = elong_sum + state_sum

            time_added = rand(Exponential(1/total_sum))
            t += time_added

            density_profile += lattice .* time_added

            RV = rand()

            if RV <= elong_sum/total_sum
                elongation_process_obc(
                    elongation_vector, internal_state_vec, posR, lattice, rates, 
                    J, k_plus, mobile, l_ribosome, track_site, jammed, tot_part,num_clust,
                    J_c, J_w, paused, c_end_pos_tracker
                )
            else
                switching_pos = change_internal_state_obc(
                    internal_state_vec, posR, lattice, rates, elongation_vector, 
                    k_minus, k_plus, mobile, paused, jammed, l_ribosome, c_end_pos_tracker, 
                    pos_first_paused_particle_bucket, time_added
                )
                return density_profile, t
            end
        end
    end
end

function gather_first_pause_distribution_profile(
    L::Int, l_ribosome::Int, track_site, init::Float64, term::Float64, elong::Float64, 
    k_minus::Float64, k_plus::Float64, n_trials::Int)
    # Preallocate arrays once:
    lattice = create_lattice_l(L, l_ribosome)             # ensure these return arrays with concrete types
    rates = create_rates_l(init, term, elong, L, l_ribosome)
    posR = positionRibosomes(lattice, l_ribosome)
    elongation_vector = get_elongation_vector_obc(lattice, rates, l_ribosome)

    mobile = [0.0]
    paused = [0.0]
    jammed = [0.0]
    internal_state_vec = get_internal_state_vec_obc(lattice, k_plus, mobile, l_ribosome)

    pos_first_paused_particle_bucket = [0.0]
    tot_part = [0.0]
    num_clust = [0.0]
    J_c = [0]
    J_w = [0]
    J = [0]
    c_end_pos_tracker = [0.0]

    density_profile = lattice
    total_time = 0
    fp_time = 0
    counter_when_passed = 0
    @inbounds for trial in 1:n_trials
        # Reset arrays for next trial efficiently:
        fill!(lattice, 0)
        fill!(elongation_vector, 0.0)
        elongation_vector[1] = rates[1]
        fill!(internal_state_vec, 0.0)
        fill!(posR, 0)
        posR[1] = 1
        mobile[1] = 0.0
        jammed[1] = 0.0
        paused[1] = 0.0

        profile, t, fp_time = run_until_first_pause_profile(
            l_ribosome, k_minus, k_plus,
            lattice, rates, posR, elongation_vector,
            mobile, jammed, paused, internal_state_vec, track_site,
            pos_first_paused_particle_bucket, tot_part, num_clust, J, J_c, J_w,
            c_end_pos_tracker 
        )
        density_profile += profile ./ t
        total_time += t
        fp_time += first_passage_time
    end

    averaged_density_profile = density_profile ./ n_trials
    averaged_total_time = total_time /n_trials
    averaged_fp_time = fp_time / n_trials 

    return averaged_density_profile, averaged_total_time,averaged_fp_time
end

function run_until_first_pause(
    l_ribosome::Int, k_minus::Float64, k_plus::Float64, deg_rate::Float64,
    lattice::Vector{Int}, rates::Vector{Float64}, posR::Vector{Int}, elongation_vector::Vector{Float64}, w_elong::Weights,
    w_state::Weights, mobile::Vector{Float64}, jammed::Vector{Float64}, paused::Vector{Float64}, internal_state_vec::Vector{Float64}, track_site::Int64,
    pos_first_paused_particle_bucket::Vector{Float64}, tot_part::Vector{Float64}, num_clust::Vector{Float64}, 
    J::Vector{Int}, J_c::Vector{Int}, J_w::Vector{Int}, c_end_pos_tracker::Vector{Float64})
    
    # Initialize time and counters.
    t = 0.0
    no_pause_before = 0    # counts the number of exit events (elongation leading to exit) before a pause
    pause_before_indicator = 0  # 1 if no exit occurs before the pause; 0 otherwise
    pos_first_paused_particle_bucket = initiate_cluster_buckets_obc(lattice)
    c_end_pos_tracker = [0.0]
    while true
        #w_elong.values .= elongation_vector # sum of all the rates in the elongation vector
        w_elong.sum    = sum(elongation_vector)
        #w_state.values .= internal_state_vec
        w_state.sum    = sum(internal_state_vec)
        total_sum = w_elong.sum  + w_state.sum + deg_rate
        time_added = rand(Exponential(1 / total_sum))
        t += time_added
        RV = rand()
        if RV <= w_elong.sum / total_sum
            # An elongation event occurs.
            elongation_process_obc(
                    elongation_vector, w_elong, internal_state_vec, posR, lattice, rates, 
                    J, k_plus, mobile, l_ribosome, track_site, jammed, tot_part,num_clust,
                    J_c, J_w, paused,c_end_pos_tracker
                    )

            # Check if a particle leaves the system.

        elseif w_elong.sum / total_sum < RV <= (w_elong.sum  + w_state.sum) / total_sum
            # An internal state change occurs.
            switching_pos = change_internal_state_obc!(
                    internal_state_vec, w_state, posR, lattice, rates, elongation_vector, 
                    k_minus, k_plus, mobile, paused, jammed, l_ribosome,c_end_pos_tracker,pos_first_paused_particle_bucket,1)

            return switching_pos,t
        else
            #reset_observables(lattice, rates, elongation_vector, internal_state_vec, posR, tot_part, mobile, jammed, paused)
            return nothing, t
        end
    end
end


function gather_first_pause_distribution(
    L::Int, l_ribosome::Int, track_site::Int64, init::Float64, term::Float64, elong::Float64, deg_rate,
    k_minus::Float64, k_plus::Float64, n_trials::Int)
    mobile = [0.0]
    paused = [0.0]
    jammed = [0.0]

    pos_first_paused_particle_bucket = [0.0]
    tot_part = [0.0]
    num_clust = [0.0]
    J_c = [0]
    J_w = [0]
    J = [0]
    c_end_pos_tracker = [0.0]

    # Preallocate arrays once:
    lattice = create_lattice_l(L,l_ribosome) # create lattice (lattice= start site(always empty) + init_region (length(l_ribosome)) + L + l_ribosome)
    rates = create_rates_l(init, term, elong, L) # rate vector
    posR = positionRibosomes(lattice) # position vector
    elongation_vector = get_elongation_vector_obc(lattice, rates) # which particles contribute to configuration
    internal_state_vec = get_internal_state_vec_obc(lattice, k_plus, mobile) # which particle contributes to configuration via its internal state
    
    w_elong  = Weights(elongation_vector)
    w_elong.values === elongation_vector
    
    w_state  = Weights(internal_state_vec)
    w_state.values === internal_state_vec

    pausing_pos = 0
    total_time = 0.0
    total_trials_time = 0
    total_trials_pos = 0
    #pause_before_total = 0   # counts the number of trials with no exit before pause

    for trial in 1:n_trials
        # Reset arrays for next trial efficiently:
        fill!(lattice, 0)
        fill!(elongation_vector, 0.0)
        elongation_vector[1] = rates[1]
        fill!(internal_state_vec, 0.0)
        fill!(posR, 0)
        posR[1] = 1
        mobile[1] = 0.0
        jammed[1] = 0.0
        paused[1] = 0.0
        J[1] = 0.0

        pos, t_trial = run_until_first_pause(
            l_ribosome, k_minus, k_plus, deg_rate,
            lattice, rates, posR, elongation_vector, w_elong, w_state,
            mobile, jammed, paused, internal_state_vec, track_site,
            pos_first_paused_particle_bucket, tot_part, num_clust, J, J_c, J_w,
            c_end_pos_tracker 
        )
        if pos !== nothing
            pausing_pos += pos
            
            total_trials_pos += 1
        end
        total_time += t_trial
        total_trials_time += 1
    end
    p_entry = total_trials_pos / n_trials
    avg_pos = (total_trials_pos > 0) ? (pausing_pos / total_trials_pos) : 0.0  
    avg_time = total_time / total_trials_time
    #probability_pause_before_exit = pause_before_total / total_trials

    return avg_pos, avg_time, p_entry
end

function run_until_end_or_hit(
    l_ribosome::Int, k_minus::Float64, k_plus::Float64,
    lattice::Vector{Int}, rates::Vector{Float64}, posR::Vector{Int}, elongation_vector::Vector{Float64},
    mobile::Vector{Float64}, jammed::Vector{Float64}, paused::Vector{Float64}, internal_state_vec::Vector{Float64}, track_site::Int64,
    pos_first_paused_particle_bucket::Vector{Float64}, tot_part::Vector{Float64}, num_clust::Vector{Float64}, 
    J::Vector{Int}, J_c::Vector{Int}, J_w::Vector{Int}, c_end_pos_tracker::Vector{Float64})::Tuple{Int64, Int64}    
    # Initialize time and counters.
    t = 0.0
    no_pause_before = 0    # counts the number of exit events (elongation leading to exit) before a pause
    pause_after_indicator = 0  # 1 if no exit occurs before the pause; 0 otherwise

    while true
        elong_sum = sum(elongation_vector)
        state_sum = sum(internal_state_vec)
        total_sum = elong_sum + state_sum
        time_added = rand(Exponential(1 / total_sum))
        t += time_added
        RV = rand()
        if RV <= elong_sum / total_sum
            # An elongation event occurs.
            elongation_process_obc(
                elongation_vector, internal_state_vec, posR, lattice, rates, 
                J, k_plus, mobile, l_ribosome, track_site, jammed, tot_part, num_clust,
                J_c, J_w, paused, c_end_pos_tracker
            )

            # Check if a particle leaves the system.
            if J[1] == 1
                no_pause_before = 1
                return no_pause_before, 0  # Mark that an exit has occurred.
            end
        else
            # A pausing event (internal state change) occurs.
            switching_pos = change_internal_state_obc(
                internal_state_vec, posR, lattice, rates, elongation_vector, 
                k_minus, k_plus, mobile, paused, jammed, l_ribosome, c_end_pos_tracker, 
                pos_first_paused_particle_bucket, time_added
            )

            # Set the indicator: if no particle left before this pause, indicator = 1.
            if no_pause_before == 0
                pause_after_indicator = 1
            end
            return 0, pause_after_indicator
            #return switching_pos, t, no_pause_before, pause_before_indicator
        end
    end
end

function gather_hit_probability(
    L::Int64, l_ribosome::Int, track_site::Int64, init::Float64, term::Float64, elong::Float64, 
    k_minus::Float64, k_plus::Float64, n_trials::Int)
    # Preallocate arrays once:
    lattice = create_lattice_l(L, l_ribosome)
    rates = create_rates_l(init, term, elong, L, l_ribosome)
    posR = positionRibosomes(lattice, l_ribosome)
    elongation_vector = get_elongation_vector_obc(lattice, rates, l_ribosome)

    mobile = [0.0]
    paused = [0.0]
    jammed = [0.0]
    internal_state_vec = get_internal_state_vec_obc(lattice, k_plus, mobile, l_ribosome)

    pos_first_paused_particle_bucket = [0.0]
    tot_part = [0.0]
    num_clust = [0.0]
    J_c = [0]
    J_w = [0]
    J = [0]
    c_end_pos_tracker = [0.0]

    pausing_pos = 0
    total_time = 0.0
    total_trials = 0
    pause_before_total = 0   # counts the number of trials with no exit before pause

    for trial in 1:n_trials
        # Reset arrays for next trial efficiently:
        fill!(lattice, 0)
        fill!(elongation_vector, 0.0)
        elongation_vector[1] = rates[1]
        fill!(internal_state_vec, 0.0)
        fill!(posR, 0)
        posR[1] = 1
        mobile[1] = 0.0
        jammed[1] = 0.0
        paused[1] = 0.0
        J[1] = 0.0

        pause_count, pause_indicator =  run_until_end_or_hit(
            l_ribosome, k_minus, k_plus,
            lattice, rates, posR, elongation_vector,
            mobile, jammed, paused, internal_state_vec, track_site,
            pos_first_paused_particle_bucket, tot_part, num_clust, 
            J, J_c, J_w, c_end_pos_tracker)

        #pausing_pos += pos
        #total_time += t_trial
        total_trials += 1
        pause_before_total += pause_count
    end

    #avg_pos = pausing_pos / total_trials
    #avg_time = total_time / total_trials
    probability_pause_before_exit = pause_before_total / total_trials

    return probability_pause_before_exit
end

function run_until_first_pause_single_particle(kp, elong, L)
    pos = 1  # Start checking from the track site

    while pos <= L + 1
        # Calculate total rate: elongation or pausing
        total_rate = elong + kp

        # Time until next event (not used here since we're tracking position)
        # Sample the next event
        if rand() < kp / total_rate
            return pos  # Pause occurred at current position
        else
            pos += 1  # Move track site forward
        end
    end

    return L + 1  # Indicate no pause within the lattice
end

function gather_first_pause_single_distribution(kp, elong, L, n_trials)
    pos_counts = zeros(Int, L + 1)  # 1 to L

    for _ in 1:n_trials
        pos = run_until_first_pause_single_particle(kp, elong, L)
        if pos <= L
            pos_counts[pos] += 1
        end
    end

    # Normalize to get probabilities
    pos_distribution = pos_counts ./ sum(pos_counts)
    return pos_distribution
end

function precompile_task()
    # Small, non-intensive parameters
    L, l_ribosome, track_site = 1, 1, 1
    deg_t, deg_g, init, term, elong = 0, 0, 0.1, 1, 1
    k₋, k₊ = 0.0001, 0.00001
    run_time, starting_t, delta_t = 1, 1, 1

    # Call your function with these dummy parameters
    Gillespie_obc(L, l_ribosome, track_site, deg_t, deg_g, init, term, elong, k₋, k₊, run_time, starting_t, delta_t; kymo=false)
end

function plot_kymo(lattice_matr, internal_state_matr, f, τ, t_vec)
    m_size = 1
    if f == 1/τ
        f *= 0.999999
    end
    plo = scatter(legend=false, xlabel = "lattice site", ylabel = "t") # initiate kymo plot
    i = 1 # track rows
    for time in t_vec # for each time step in time vector
      lattice = lattice_matr[i, 1:end] # choose the lattice corresponding to the time
      internal_state_vec = internal_state_matr[i, 1:end]# choose the lattice corresponding to the internal state
      t = time # set y axis
      for index in eachindex(lattice[1:end-1]) # go through lattice and assign points for specifc particle classes
        # mobile 
        if lattice[index] == 1 && lattice[index+1] == 0 && internal_state_vec[index] == f #green square
          scatter!([index],[t], markercolor = "lightgreen", markersize = m_size, markershape = :rect, markerstrokewidth=0)
        end
  
        # jammed by mobile or paused
        if lattice[index] == 1 && lattice[index+1] == 1 && internal_state_vec[index] == f
          scatter!([index], [t], markercolor="lightblue", markersize = m_size, markershape = :rect, markerstrokewidth=0)
        end
  
        if lattice[index] == 1 && internal_state_vec[index] == 1/τ
          scatter!([index], [t], markercolor="black", markersize = m_size, markershape = :rect, markerstrokewidth=0)
        end
        
        # add holes as color
        # if lattice[index] == 0
        #   scatter!([index], [t], markercolor="lightgrey", markersize = 2.5, markershape = :rect, markerstrokewidth=0)
        # end

        # periodic boundary conditions
      end

      if lattice[end] == 1 && lattice[1] == 1 && internal_state_vec[end] == f
        scatter!([length(lattice)], [t], markercolor="lightgrey", markersize = m_size, markershape = :rect, markerstrokewidth=0)
      end

      if lattice[end] == 1 && lattice[1] == 0 && internal_state_vec[end] == f #green square
        scatter!([length(lattice)],[t], markercolor = "lightgreen", markersize = m_size, markershape = :rect, markerstrokewidth=0)
      end

      if lattice[end] == 1 && internal_state_vec[end] == 1/τ
        scatter!([length(lattice)], [t], markercolor="black", markersize = m_size, markershape = :rect, markerstrokewidth=0)
      end

      i+=1

    end
    
    return plo
end


function plot_kymo2(lattice_matr, internal_state_matr, f, τ, t_vec, l_ribosome, track_site)
    m_size = 2
    if f == 1/τ
        f *= 0.999999
    end
    plo = scatter(legend=false, size =(2000, 1200), xaxis =false, yaxis = false, xlabel=false, ylabel=false) # initiate kymo plot
    num_sites = size(lattice_matr, 2)
    num_times = length(t_vec)

    x_vals = Int[]
    y_vals = Float64[]
    colors = String[]

    for i in 1:num_times
        lattice = lattice_matr[i, :]
        internal_state_vec = internal_state_matr[i, :]
        t = t_vec[i]
        for index in 1:num_sites
            if lattice[index] == 1
                # Adjusted calculation for particle_sites
                start_site = max(1, index - (track_site - 1))
                end_site = min(num_sites, index + (l_ribosome - track_site))
                particle_sites = start_site:end_site

                # Determine the type of particle for coloring
                if internal_state_vec[index] == f
                    marker_color = "lightgreen"
                elseif internal_state_vec[index] == 1/τ
                    marker_color = "black"
                else
                    marker_color = "lightblue"
                end

                # Collect data for plotting
                for site in particle_sites
                    push!(x_vals, site)
                    push!(y_vals, t)
                    push!(colors, marker_color)
                end
            end
        end
    end

    # Plot all points at once
    scatter!(x_vals, y_vals, markercolor=colors, markersize=m_size, markershape=:rect, markerstrokewidth=0)

    return plo
end

function plot_kymo2(lattice_matr, internal_state_matr, f, τ, t_vec, l_ribosome, track_site;
                    fig_size::Tuple{Int,Int} = (800, 600),
                    fontfamily::AbstractString = "Cambria Math",
                    guidefontsize::Int = 18,
                    tickfontsize::Int = 12,
                    legendfontsize::Int = 12)

    m_size = 2.5
    if f == 1/τ
        f *= 0.999999
    end

    # base plot with fonts
    plo = scatter(
        legend=false,
        size=fig_size,
        fontfamily=fontfamily,                # global family
        guidefont=font(guidefontsize, fontfamily),   # axes labels
        tickfont=font(tickfontsize, fontfamily),     # tick labels
        legendfont=font(legendfontsize, fontfamily)  # legend (if used)
    )

    num_sites = size(lattice_matr, 2)
    num_times = length(t_vec)

    x_vals = Int[]
    y_vals = Float64[]
    colors = String[]

    for i in 1:num_times
        lattice = lattice_matr[i, :]
        internal_state_vec = internal_state_matr[i, :]
        t = t_vec[i]
        for index in 1:num_sites
            if lattice[index] == 1
                start_site = max(1, index - (track_site - 1))
                end_site   = min(num_sites, index + (l_ribosome - track_site))
                marker_color =
                    internal_state_vec[index] == f    ? "lightgreen" :
                    internal_state_vec[index] == 1/τ  ? "black"      :
                                                         "lightblue"
                for site in start_site:end_site
                    push!(x_vals, site)
                    push!(y_vals, t)
                    push!(colors, marker_color)
                end
            end
        end
    end

    scatter!(
        x_vals, y_vals;
        markercolor=colors,
        markersize=m_size,
        markershape=:rect,
        markerstrokewidth=0
    )

    return plo
end

using CairoMakie, Colors

function plot_kymo2_heatmap(lattice_matr, internal_state_matr, f, τ, t_vec,
                            l_ribosome::Int, track_site::Int;
                            fig_size::Tuple{Int,Int} = (1200, 600),
                            yflip::Bool = true,
                            atol::Real = 0.0,
                            rasterize_body::Bool = true)

    # --- normalize orientation (rows=time) ---
    t = Float64.(t_vec)
    A = Array(lattice_matr)
    S = Array(internal_state_matr)

    if size(A, 1) != length(t) && size(A, 2) == length(t)
        A = A'
        S = S'
    end

    # common off-by-one fix (e.g. extra t0 frame)
    if size(A, 1) == length(t) + 1
        A = A[2:end, :]; S = S[2:end, :]
    elseif size(A, 1) == length(t) - 1
        t = t[1:end-1]
    end

    @assert size(A,1) == length(t) "time length must match number of rows in lattice_matr."
    @assert size(A) == size(S)     "lattice_matr and internal_state_matr must be same size."

    T, L = size(A)

    # --- build class matrix K: 0=bg, 1=unpaused, 2=paused, 3=other ---
    K = zeros(Int, T, L)
    f_unpause = float(f)
    f_pause   = 1 / float(τ)

    @inbounds for i in 1:T, j in 1:L
        if A[i, j] == 1
            s = max(1, j - (track_site - 1))
            e = min(L, j + (l_ribosome - track_site))
            v = S[i, j]
            code = (atol == 0 && v == f_unpause) ? 1 :
                   (atol == 0 && v == f_pause)   ? 2 :
                   (atol > 0  && isapprox(v, f_unpause; atol=atol, rtol=0)) ? 1 :
                   (atol > 0  && isapprox(v, f_pause;   atol=atol, rtol=0)) ? 2 : 3
            K[i, s:e] .= code
        end
    end

    # --- colors: consistent RGBA{Float64} types ---
    cmap = RGBA{Float64}.([
        RGBA(1, 1, 1, 0.0),      # 0 background (transparent)
        colorant"lightgreen",    # 1 unpaused
        colorant"black",         # 2 paused
        colorant"lightblue"      # 3 other
    ])

    # --- figure/axis ---
    fig = Figure(size = fig_size)
    ax  = Axis(fig[1, 1]; xlabel = "Lattice Site", ylabel = "Time", yreversed = yflip)

    # --- build edges (required: length(x) = L+1, length(y) = T+1) ---
    xedges = collect(0.5:1:(L + 0.5))                    # L+1
    ymids  = 0.5 .* (t[1:end-1] .+ t[2:end])             # T-1
    y0     = t[1]  - (t[2]  - t[1])  / 2
    yend   = t[end] + (t[end] - t[end-1]) / 2
    yedges = vcat(y0, ymids, yend)                       # T+1

    # --- heatmap with edges; transpose K so edges match (ni=L, nj=T) ---
    Makie.heatmap!(ax, xedges, yedges, Float32.(K)';     # <-- K'
                    colormap = cmap,
                    colorrange = (0, 3),
                    interpolate = false,
                    rasterize = rasterize_body)

    return fig
end


function cluster_len_num(time_added, cluster_num_bucket, cluster_len_bucket, lattice, num_clust)
    holes = findall(x->x==0 , lattice) # findall the position of all holes on the lattice
    p_clusters = [holes[item+1] - holes[item] for item in eachindex(holes[1:end-1])] .- 1 # distance between holes
    filter!(x-> x >= 1, p_clusters) # disregard all neighboring holes
    for item in p_clusters # add time to specific cluster length
        cluster_len_bucket[item] +=  time_added
    end
    cluster_num_bucket[num_clust[1]] += time_added
end

function check_num_cluster(lattice, num_clust)
    test = true
    holes = findall(x->x==0 , lattice) # findall the position of all holes on the lattice
    p_clusters = [holes[item+1] - holes[item] for item in eachindex(holes[1:end-1])] .- 1
    filter!(x-> x >= 1, p_clusters)
    num_clust2 = length(p_clusters)

    if num_clust2 != num_clust[1]
        test=false
    end
    return test 
end

function number_of_cluster(pos_paused_particles, lattice, paused)
    num_cluster = paused[1] # set number of clusters
    # pushfirst!(pos_paused_particles, findlast(x->x==0, lattice[1:pos_paused_particles[1]])) #add a zero to paused particles so that first cluster is also counted
    for i in eachindex(pos_paused_particles[1:end-1])
        if sum(lattice[pos_paused_particles[i]+1:pos_paused_particles[i+1]]) == pos_paused_particles[i+1] - pos_paused_particles[i]
            num_cluster -= 1
        end
    end
    return num_cluster
end

function paused_particle_distribution_single_cluster(
    lattice, internal_state_vec, paused_dist_cluster, 
    pos_paused_particles, k₊, time_added
    )
    
    if isempty(pos_paused_particles) == false # check if there are any paused particles
        for item in reverse(pos_paused_particles) 
            count_backwards(lattice, internal_state_vec, paused_dist_cluster, item, k₊, time_added)
        end
    else
        paused_dist_cluster[-1] += time_added
    end
end

function count_backwards(lattice, internal_state_vec, paused_dist_cluster, start_pos, k₊, time_added)
    count = 0
    stop = false
    while stop == false
        if lattice[start_pos-1-count] == 1 && internal_state_vec[start_pos-1-count] == k₊
            count += 1
        else
            paused_dist_cluster[count] += time_added
            stop = true 
        end
    end
    return count
end

function sort_clusters_length_and_dist(
    lattice, paused, internal_state_vec, paused_dist_cluster, 
    k₊, time_added, k₋)

    pos_paused_particles = findall(x->x==k₋, internal_state_vec)
    # num_clusters = number_of_cluster(pos_paused_particles, lattice, paused)
    paused_particle_distribution_single_cluster(lattice, internal_state_vec, paused_dist_cluster, pos_paused_particles, k₊, time_added)
end

# cluster_buckets = initiate_cluster_buckets_obc(L, l_ribosome)
function sort_clusters_obc(cluster_buckets, lattice, l_ribosome)
    particles = findall(x->x==1, lattice[1:end-l_ribosome]) # findall the position of all holes on the lattice
    p_clusters = [particles[item+1] - particles[item] for item in eachindex(particles[1:end-1])]
    cluster_counter = 1
    for item in p_clusters
        if item == l_ribosome
            cluster_counter += 1
        else

            if cluster_counter != 1 # just single particles
                cluster_buckets[cluster_counter] += 1
                cluster_counter = 1
            else
                cluster_buckets[1] += 1
            end
    
        end
    end

    if lattice[particles[1]] == 1 && lattice[particles[1]+1] == 0 
        cluster_buckets[1] += 1
    end
end

function cluster(lattice)
    holes = findall(x->x==0 , lattice) # findall the position of all holes on the lattice
    particles = findall(x->x ==1, lattice) # findall the position of particles
    #count the number of sites between the holes and the particles
    p_clusters = [holes[item+1] - holes[item] for item in eachindex(holes[1:end-1])] .- 1
    h_clusters = [particles[item+1] - particles[item] for item in eachindex(particles[1:end-1])] .-1
    # pbc means that cluster can be extended after L+1
    pbc_p_cluster = length(lattice) - holes[end] + holes[1] - 1
    pbc_h_cluster = length(lattice) - particles[end] + particles[1] -1 
    push!(p_clusters, pbc_p_cluster)
    push!(h_clusters, pbc_h_cluster) 
    # a cluster must be more then one particle
    filter!(x-> x >= 1, p_clusters)
    filter!(x-> x >= 1, h_clusters)
    # error catch : no cluster on the lattice returns 0 instead of empty vector
    if isempty(p_clusters) == true
        p_clusters = [0.0]
    elseif isempty(h_clusters) == true
        h_clusters = [0.0]
    end

    return h_clusters, p_clusters
end

function sort_clusters(h_clusters, h_cluster_buckets, p_clusters, p_cluster_buckets, time_added, lattice)
    test = false
    h_clusters_classes = sort(unique(h_clusters))
    p_clusters_classes = sort(unique(p_clusters))

    number_of_h_clusters = [count(x->x==item, h_clusters) for item in h_clusters_classes] .* time_added
    number_of_p_clusters = [count(x->x==item, p_clusters) for item in p_clusters_classes] .* time_added

    if sum(number_of_p_clusters .* p_clusters_classes) != sum(lattice)
        test = true
    end

    for cluster in eachindex(h_clusters_classes)
        h_cluster_buckets[h_clusters_classes[cluster]] += number_of_h_clusters[cluster]
    end

    for cluster in eachindex(p_clusters_classes)
        p_cluster_buckets[p_clusters_classes[cluster]] += number_of_p_clusters[cluster]
    end
    
    return test
end

function initiate_cluster_buckets_obc(lattice)

    # h_cluster_buckets = Dict{Float64, Float64}()
    p_cluster_buckets = Dict{Float64, Float64}()

    # for j in 1:20
    #     p_cluster_buckets[j] = Dict()
    #     for i in 0:L/l_ribosome
    #         p_cluster_buckets[j][i] = 0.0
    #     end
    # end
    for i in 0:length(lattice)
        p_cluster_buckets[i] = 0.0
    end
    p_cluster_buckets[-1.0] = 0.0 
    # holes = 1-ρ
    # for i in 1:holes*L
        # h_cluster_buckets[i] = 0
    # end

    return p_cluster_buckets
end

function clean_clusters(clusters)

    for key in keys(clusters)
        if sum(values(clusters[key])) == 0
            delete!(clusters, key)
        end
    end

    # for key in keys(clusters)
    #     filter!(p->p.second != 0, clusters[key])
    # end

end

function normalize_cluster(cluster_data)
    total_time_α = []
    time_in_cluster = []

    for dict in cluster_data
        dic = Dict()
        
        total_time = 0
        for key in keys(dict) 
            total_time += sum(values(dict[key])) 
        end
        for key in keys(dict)
            dic[key] = sum(values(dict[key]))/total_time
        end

        push!(time_in_cluster , dic)
        push!(total_time_α, total_time)
        # push!(norm_dist_clust, dic2)
    end
    test = cluster_data

    return total_time_α, time_in_cluster
end


function Gillespie_obc_test(L, l_ribosome, track_site, deg_t, deg_g, init, term, elong, k₋, k₊, run_time, starting_t, delta_t; kymo=false)
    
    if k₋ == k₊
        k₊ *= 0.99999999
    end

    # For consistency with your existing code, we still define number_of_measurements.
    number_of_measurements = Int64((run_time - starting_t) / delta_t)

    # Observable counters and vectors
    mobile = [0.0] 
    paused = [0.0]
    jammed = [0.0]
    tot_part = [0.0]
    tot_weighted = [0.0]
    tot_weighted_p = [0.0]
    tot_weighted_unp = [0.0]
    tot_vec_p = Vector{Float64}(undef, number_of_measurements)
    tot_vec_unp = Vector{Float64}(undef, number_of_measurements)
    tot_vec = Vector{Float64}(undef, number_of_measurements)
    num_clust = [0.0]
    unpausing_counter = [0.0]
    unpausing_rate  = [0.0]
    pausing_counter = [0.0]
    pausing_rate_1 = [0.0]
    t_c_vec = Vector{Float64}(undef, number_of_measurements)
    t_w_vec = Vector{Float64}(undef, number_of_measurements)  
    time_vec = Vector{Float64}(undef, number_of_measurements)
    mobile_vec = Vector{Float64}(undef, number_of_measurements)
    mobile_weighted = [0.0]
    jammed_vec = Vector{Float64}(undef, number_of_measurements)
    jammed_weighted = [0.0]
    paused_vec = Vector{Float64}(undef, number_of_measurements)
    paused_weighted = [0.0]
    current = Vector{Float64}(undef, number_of_measurements)
    w_current = Vector{Float64}(undef, number_of_measurements)
    c_current = Vector{Float64}(undef, number_of_measurements)
    J_w = [0] # TASEP current weighted by time
    J_c = [0] # pTASEP current
    J = [0]   # current counter (e.g., particles leaving last site)
    J_c_end = [0]
    c_end_pos_tracker = [0.0]

    lattice = create_lattice_l(L, l_ribosome)
    rates = create_rates_l(init, term, elong, L, l_ribosome)
    posR = positionRibosomes(lattice, l_ribosome)
    elongation_vector = get_elongation_vector_obc(lattice, rates, l_ribosome)
    internal_state_vec = get_internal_state_vec_obc(lattice, k₊, mobile, l_ribosome)
    cluster_num_bucket_p = initiate_cluster_buckets_obc(lattice)
    cluster_num_bucket_unp = initiate_cluster_buckets_obc(lattice)
    cluster_len_bucket = initiate_cluster_buckets_obc(lattice)
    paused_dist_bucket = initiate_cluster_buckets_obc(lattice)
    pos_first_paused_particle_bucket = initiate_cluster_buckets_obc(lattice)

    # ----- Automatic Steady-State Detection -----
    # We use sliding windows of duration delta_t and update an EMA for both the current and the density.
    steady_state_detected = false
    window_time = 0.0          # time accumulator for the current window
    J_window = 0.0             # accumulated current in the window (we record the increment)
    tot_weighted_window = 0.0  # accumulated weighted density in the window
    window_count = 0
    consecutive_windows = 0    # count how many consecutive windows satisfy the criterion
    alpha = 0.05                # EMA factor (adjust for sensitivity)
    threshold = 0.025           # 1% relative difference threshold
    steady_window_required = 3 # require 3 consecutive windows meeting the criterion
    eps = 1e-10                # small number to avoid division by zero
    ema_current = 0.0
    ema_density = 0.0
    t = 0.0  # simulation time
    prev_J = J[1]  # to capture increments in current

    println("Starting automatic steady-state detection...")
    while !steady_state_detected && t < run_time
        # Compute event rates
        elong_sum = sum(elongation_vector)
        state_sum = sum(internal_state_vec)
        if paused[1] == 0
            deg_rate = deg_g
        else
            deg_rate = deg_t
        end
        total_sum = elong_sum + state_sum + deg_rate
        RV = rand()
        time_added = rand(Exponential(1/total_sum))
        
        # Choose event based on RV and update the state
        if RV <= elong_sum/total_sum
            elongation_process_obc(
                elongation_vector, internal_state_vec, posR, lattice, rates, 
                J, k₊, mobile, l_ribosome, track_site, jammed, tot_part, num_clust,
                J_c, J_w, paused, c_end_pos_tracker
            )
        elseif RV <= (elong_sum + state_sum) / total_sum
            change_internal_state_obc(
                internal_state_vec, posR, lattice, rates, elongation_vector, 
                k₋, k₊, mobile, paused, jammed, l_ribosome, c_end_pos_tracker, pos_first_paused_particle_bucket, delta_t
            )
        else
            #reset_observables(lattice, rates, elongation_vector, internal_state_vec, posR, tot_part, mobile, jammed, paused)
        end
        
        # Compute the increment in the current from this event and accumulate
        delta_J = J[1] - prev_J
        J_window += delta_J
        prev_J = J[1]
        
        # Accumulate density (weighted by the event time)
        tot_weighted_window += tot_part[1] / length(posR) * time_added
        
        window_time += time_added
        t += time_added
        
        # When the current window is complete, compute the averages and update the EMA.
        if window_time >= 10* minimum([1/k₋, 1/k₊, 1/init])
            avg_current = J_window / window_time
            avg_density = tot_weighted_window / window_time

            if window_count == 0
                ema_current = avg_current
                ema_density = avg_density
            else
                ema_current = alpha * avg_current + (1 - alpha) * ema_current
                ema_density = alpha * avg_density + (1 - alpha) * ema_density
            end
            
            diff_current = abs(avg_current - ema_current) / (ema_current + eps)
            diff_density = abs(avg_density - ema_density) / (ema_density + eps)
            
            #println("Window ", window_count, ": avg_current = ", avg_current, ", ema_current = ", ema_current,
            #        ", diff_current = ", diff_current)
            #println("Window ", window_count, ": avg_density = ", avg_density, ", ema_density = ", ema_density,
            #        ", diff_density = ", diff_density)
            
            if diff_current < threshold && diff_density < threshold
                consecutive_windows += 1
            else
                consecutive_windows = 0
            end
            
            if consecutive_windows >= steady_window_required
                steady_state_detected = true
                println("Steady state detected at simulation time t = ", t)
            end
            
            # Reset window accumulators for the next window
            window_time = 0.0
            J_window = 0.0
            tot_weighted_window = 0.0
            window_count += 1
        end
    end
    # After steady state is detected, reset counters before starting measurements
    J[1] = 0.0
    J_w[1] = 0.0
    J_c[1] = 0.0
    deg_counter_paused = [0.0]
    deg_counter_non_paused = [0.0]
    pausing_counter[1] = 0.0
    pausing_rate_1[1] = 0.0
    unpausing_counter[1] = 0.0
    unpausing_rate[1] = 0.0
    c_end_pos_tracker[1] = 0.0
    J_c_end[1] = 0.0

    t_c = 0.0
    t_w = 0.0
    t = 0.0  # reset measurement time

    # ----- Measurement Phase -----
    for i in 1:number_of_measurements
        while t <= delta_t
            elong_sum = sum(elongation_vector)
            state_sum = sum(internal_state_vec)
            
            if paused[1] == 0
                deg_rate = deg_g
            else
                deg_rate = deg_t
            end
            
            total_sum = elong_sum + state_sum + deg_rate
            RV = rand()
            time_added = rand(Exponential(1/total_sum))
            
            if paused[1] == 0 && c_end_pos_tracker[1] == 0
                t_w += time_added
                tot_weighted_unp[1] += tot_part[1] / length(posR) * time_added
            elseif paused[1] >= 1 || c_end_pos_tracker[1] != 0
                t_c += time_added
                tot_weighted_p[1] += tot_part[1] / length(posR) * time_added
            end
            
            t += time_added
            mobile_weighted[1] += mobile[1] / length(posR) * time_added
            paused_weighted[1] += paused[1] / length(posR) * time_added
            jammed_weighted[1] += jammed[1] / length(posR) * time_added
            tot_weighted[1] += tot_part[1] / length(posR) * time_added
            
            find_all_cluster_lengths_obc(
                lattice, 
                cluster_len_bucket, 
                cluster_num_bucket_unp,
                cluster_num_bucket_p, 
                time_added, 
                num_clust,
                paused,
                l_ribosome,
                c_end_pos_tracker
            )
            
            paused_dist_bucket[paused[1]] += time_added
            
            if RV <= elong_sum/total_sum
                elongation_process_obc(
                    elongation_vector, internal_state_vec, posR, lattice, rates, 
                    J, k₊, mobile, l_ribosome, track_site, jammed, tot_part, num_clust,
                    J_c, J_w, paused, c_end_pos_tracker
                )
            elseif RV <= (elong_sum + state_sum) / total_sum
                change_internal_state_obc(
                    internal_state_vec, posR, lattice, rates, elongation_vector, 
                    k₋, k₊, mobile, paused, jammed, l_ribosome, c_end_pos_tracker, pos_first_paused_particle_bucket, delta_t
                )
            else
                if paused[1] == 0 
                    deg_counter_non_paused[1] += 1
                else
                    deg_counter_paused[1] += 1
                end
                reset_observables(lattice, rates, elongation_vector, internal_state_vec, posR, tot_part, mobile, jammed, paused)
            end
        end
        
        current[i] = J[1] / t
        w_current[i] = J_w[1] / t
        c_current[i] = J_c[1] / t
        mobile_vec[i] = mobile_weighted[1] / t
        paused_vec[i] = paused_weighted[1] / t
        jammed_vec[i] = jammed_weighted[1] / t
        tot_vec[i] = tot_weighted[1] / t
        tot_vec_unp[i] = tot_weighted_unp[1] / t
        tot_vec_p[i] = tot_weighted_p[1] / t
        time_vec[i] = t
        t_w_vec[i] = t_w
        t_c_vec[i] = t_c
        
        t = 0.0
        t_w = 0.0
        t_c = 0.0
        J[1] = 0.0
        J_w[1] = 0.0
        J_c[1] = 0.0
        mobile_weighted[1] = 0.0
        paused_weighted[1] = 0.0
        jammed_weighted[1] = 0.0
        tot_weighted[1] = 0.0
    end

    J_mean::Float64 = mean(current)
    J_w_mean::Float64 = mean(w_current)
    J_c_mean::Float64 = mean(c_current)

    ρ_vec::Vector{Float64} = vcat(
        mean(mobile_vec), 
        mean(paused_vec), 
        mean(jammed_vec), 
        mean(tot_vec),  
        mean(tot_vec_unp),
        mean(tot_vec_p)
    )

    time_in_w = mean(t_w_vec)
    time_in_p = mean(t_c_vec)
    total_time = mean(time_vec)

    paused_particle_distribution = Dict(key => value/total_time for (key, value) in paused_dist_bucket)
    cluster_length_distribution = Dict(key => value/total_time for (key, value) in cluster_len_bucket)
    cluster_number_distribution_unp = Dict(key => value/total_time for (key, value) in cluster_num_bucket_unp)
    cluster_number_distribution_p = Dict(key => value/total_time for (key, value) in cluster_num_bucket_p)

    return [J_w_mean, J_c_mean, J_mean], ρ_vec, paused_particle_distribution, cluster_length_distribution, cluster_number_distribution_unp, cluster_number_distribution_p, pos_first_paused_particle_bucket, [time_in_w, time_in_p, total_time]
end

#=
begin
    L, l_ribosome, track_site = 10, 2 ,1
    ell = l_ribosome
    ϵ = 1.0
    init, term, elong = 0.5 , ϵ, ϵ
    deg_t = deg_g = 0
    number_of_measurements = 1
    k₋ = 0.5#/ ϵ
    k₊ = 0.1 #/ ϵ
    # number of different particle classes
    mobile = [0.0] 
    paused = [0.0]
    jammed = [0.0]

    tot_part = [0.0]
    tot_weighted = [0.0]
    tot_weighted_p = [0.0]
    tot_weighted_unp = [0.0]

    tot_vec_p = Vector{Float64}(undef,number_of_measurements)
    tot_vec_unp = Vector{Float64}(undef,number_of_measurements)
    tot_vec = Vector{Float64}(undef,number_of_measurements)

    num_clust = [0.0]

    unpausing_counter = [0.0]
    unpausing_rate  = [0.0]

    pausing_counter = [0.0]
    pausing_rate_1 = [0.0]

    t_c_vec = Vector{Float64}(undef, number_of_measurements)
    t_w_vec = Vector{Float64}(undef, number_of_measurements)  
    time_vec = Vector{Float64}(undef, number_of_measurements)

    # ρ_vec = Vector{Float64}(undef,number_of_measurements) # density vector
    # ρ_weighted = [0.0] # density after weighted by time step the configuration

    mobile_vec = Vector{Float64}(undef, number_of_measurements) # mobile vector
    mobile_weighted = [0.0] # mobile fraction of particles

    jammed_vec= Vector{Float64}(undef,number_of_measurements)# jammed vector
    jammed_weighted = [0.0]# jammed fraction of particles

    paused_vec = Vector{Float64}(undef,number_of_measurements)# paused vector
    paused_weighted = [0.0]# paused fraction of particles

    current = Vector{Float64}(undef,number_of_measurements) # current vector
    w_current = Vector{Float64}(undef,number_of_measurements)
    c_current = Vector{Float64}(undef,number_of_measurements)

    J_w = [0] # TASEP current
    J_c = [0] # pTASEP current
    J = [0] # current counter at last site
    J_c_end = [0]

    zero_to_one = [0]
    one_to_zero = [0]


    lattice = create_lattice_l(L, l_ribosome) # create lattice (lattice= start site(always empty) + init_region (length(l_ribosome)) + L + l_ribosome)
    rates = create_rates_l(init, term, elong, L) # rate vector
    posR = positionRibosomes(lattice) # position vector
    elongation_vector = get_elongation_vector_obc(lattice, rates) # which particles contribute to configuration
    internal_state_vec = get_internal_state_vec_obc(lattice, k₊, mobile) # which particle contributes to configuration via its internal state
    c_end_pos_tracker = [0.0]
    w_elong  = Weights(elongation_vector)
    w_elong.values === elongation_vector
    
    w_state  = Weights(internal_state_vec)
    w_state.values === internal_state_vec
    
    cluster_num_bucket_p = initiate_cluster_buckets_obc(lattice)
    cluster_num_bucket_unp = initiate_cluster_buckets_obc(lattice)
    cluster_len_bucket = initiate_cluster_buckets_obc(lattice)
    paused_dist_bucket = initiate_cluster_buckets_obc(lattice)
    pos_first_paused_particle_bucket = initiate_cluster_buckets_obc(lattice)

    t = 0.0 # reset time counter
end

tmp = []
tmp_tot = []

t_c = 0.0
t_w = 0.0
lattice
t = 0.0

begin
     # Compute event rates
    w_elong.sum = sum(elongation_vector)

    #w_state.values .= internal_state_vec
    w_state.sum    = sum(internal_state_vec)

    if paused[1] == 0
        deg_rate = deg_g
    else
        deg_rate = deg_t
    end

    total_sum = w_elong.sum + w_state.sum  + deg_rate
    RV = rand()
    time_added = rand(Exponential(1/total_sum))
    t += time_added
    
    # Choose event based on RV and update the state
    if RV <= w_elong.sum /total_sum
        elongation_process_obc(
            elongation_vector, w_elong, internal_state_vec, posR, lattice, rates, 
            J, k₊, mobile, l_ribosome, track_site, jammed, tot_part, num_clust,
            J_c, J_w, paused, c_end_pos_tracker
        )
    elseif w_elong.sum /total_sum < RV <= (w_elong.sum  + w_state.sum) / total_sum 
        change_internal_state_obc!(
            internal_state_vec, w_state, posR, lattice, rates, elongation_vector, 
            k₋, k₊, mobile, paused, jammed, l_ribosome, c_end_pos_tracker, pos_first_paused_particle_bucket, time_added
        )
    else
        #reset_observables(lattice, rates, elongation_vector, internal_state_vec, posR, tot_part, mobile, jammed, paused)
    end
    print("total particles =  $(tot_part[1])\n")
    print("mobile =  $(mobile[1])\n")
    print("paused =  $(paused[1])\n")
    print("jammed =  $(jammed[1])\n")
    print("J =  $(J[1])\n")
    print("J_w =  $(J_w[1])\n")
    print("J_c =  $(J_c[1])\n")
    hcat(posR, internal_state_vec, elongation_vector)

end

lattice


plot(posR)

J

for i in 1:500
    
    begin
        L, l_ribosome, track_site= 1000, 10, 1  
        init, elong, term = 0.02, 1, 10
        k₊, k₋ = 0.0, 0.01
        number_of_measurements = 1
        # number of different particle classes
        mobile = [0.0] 
        paused = [0.0]
        jammed = [0.0]

        tot_part = [0.0]
        tot_weighted = [0.0]
        tot_weighted_p = [0.0]
        tot_weighted_unp = [0.0]

        tot_vec_p = Vector{Float64}(undef,number_of_measurements)
        tot_vec_unp = Vector{Float64}(undef,number_of_measurements)
        tot_vec = Vector{Float64}(undef,number_of_measurements)

        num_clust = [0.0]

        unpausing_counter = [0.0]
        unpausing_rate  = [0.0]

        pausing_counter = [0.0]
        pausing_rate_1 = [0.0]

        t_c_vec = Vector{Float64}(undef, number_of_measurements)
        t_w_vec = Vector{Float64}(undef, number_of_measurements)  
        time_vec = Vector{Float64}(undef, number_of_measurements)

        # ρ_vec = Vector{Float64}(undef,number_of_measurements) # density vector
        # ρ_weighted = [0.0] # density after weighted by time step the configuration

        mobile_vec = Vector{Float64}(undef, number_of_measurements) # mobile vector
        mobile_weighted = [0.0] # mobile fraction of particles

        jammed_vec= Vector{Float64}(undef,number_of_measurements)# jammed vector
        jammed_weighted = [0.0]# jammed fraction of particles

        paused_vec = Vector{Float64}(undef,number_of_measurements)# paused vector
        paused_weighted = [0.0]# paused fraction of particles

        current = Vector{Float64}(undef,number_of_measurements) # current vector
        w_current = Vector{Float64}(undef,number_of_measurements)
        c_current = Vector{Float64}(undef,number_of_measurements)

        J_w = [0] # TASEP current
        J_c = [0] # pTASEP current
        J = [0] # current counter at last site
        J_c_end = [0]

        zero_to_one = [0]
        one_to_zero = [0]


        lattice = create_lattice_l(L, l_ribosome) # create lattice (lattice= start site(always empty) + init_region (length(l_ribosome)) + L + l_ribosome)
        rates = create_rates_l(init, term, elong, L) # rate vector
        posR = positionRibosomes(lattice) # position vector
        elongation_vector = get_elongation_vector_obc(lattice, rates) # which particles contribute to configuration
        internal_state_vec = get_internal_state_vec_obc(lattice, k₊, mobile) # which particle contributes to configuration via its internal state
        c_end_pos_tracker = [0.0]

         w_elong  = Weights(elongation_vector)
        w_elong.values === elongation_vector
    
        w_state  = Weights(internal_state_vec)
        w_state.values === internal_state_vec

        cluster_num_bucket_p = initiate_cluster_buckets_obc(lattice)
        cluster_num_bucket_unp = initiate_cluster_buckets_obc(lattice)
        cluster_len_bucket = initiate_cluster_buckets_obc(lattice)
        paused_dist_bucket = initiate_cluster_buckets_obc(lattice)
        pos_first_paused_particle = initiate_cluster_buckets_obc(lattice)

        t = 0.0 # reset time counter
    end
    while lattice[L] == 0 
        begin
            elong_sum, state_sum, deg_rate = sum(elongation_vector), 0, 0
            total_sum = elong_sum + state_sum + deg_rate
            time_added = rand(Exponential(1/total_sum))
            t += time_added
            elongation_process_obc(
                elongation_vector, w_elong, internal_state_vec, posR, lattice, rates, 
                J, k₊, mobile, l_ribosome, track_site, jammed, tot_part,num_clust,
                J_c, J_w, paused,c_end_pos_tracker
                )
            #tmp = hcat(posR, internal_state_vec)
            #if c_end_pos_tracker[1] != 0.0
            #    #break
            #end
            #push!(tmp, hcat(posR, internal_state_vec))
            #push!(tmp_cend, c_end_pos_tracker[1])
            #println(c_end_pos_tracker[1])
            #println(J_c)
            #println(c_end_pos_tracker)
            #println(J_w)
            hcat(posR, internal_state_vec)
            #pos_first_paused_particle
        end
    end
    push!(tmp, t)
    push!(tmp_tot, tot_part[1])
end

(L)/(1-init) / mean(tmp)

bar(tmp)
lattice

using Plots
bar(lattice)
tot_part

hcat(posR, internal_state_vec)
J_c  + J_w
J

print("hello")
=#