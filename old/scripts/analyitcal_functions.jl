using Pkg
using LsqFit
using Statistics

begin 
    # γ_const = γ * ρ * Φ_const * (1-Φ_const)
    #=
    Parameters:
        α = initiation rate
        β = termination rate
        γ = hopping rate
        k₊ = rate of entering paused state
        k₋ = rate at which particles leave
        K_A = rate of association
    =#
    P0(k₋, k₊, ρ,L) = fa(k₋, k₊)^(ρ*L)
    fa(k₋, k₊) =  k₋ / (k₋+k₊)

    fa_hat(k₋, k₊,ρ, ϵ) = (k₋/(k₋ + k₊ + ϵ *ρ*fp(k₋, k₊)) )
    fa_hat_ext(k₋, k₊,ρ, ϵ, l_ribosome) = (k₋/(k₋ + k₊ + ϵ *ρ*1 *fp(k₋, k₊)) )
    fp(k₋, k₊) = k₊ /(k₋+k₊)

    α_crit_hat_fct(k₋, k₊, ϵ) = ρ_max_fct(k₋,k₊, ϵ) * ϵ * fa_hat(k₋, k₊, ρ_max_fct(k₋,k₊, ϵ), ϵ)
    α_crit_hat_ext(k₋, k₊, ϵ, ρ_num) = ρ_num * ϵ * fa_hat(k₋, k₊, ρ_num, ϵ) / (l_ribosome - ρ_num*l_ribosome + ρ_num)
    
    ρ_LD_Wang_ext(init,elong,l_ribosome,kp,km) = init * l_ribosome /(elong*fa_hat_ext(km, kp,ρ, ϵ, l_ribosome) + init*(l_ribosome+1))    

    ρ_max_fct(k₋,k₊, γ) = -(k₋+k₊)/(γ*fp(k₋, k₊)) + sqrt( ((k₋+k₊) / (γ*fp(k₋,k₊)))^2 + (k₋+k₊)/(γ*fp(k₋, k₊))) 
    
    Wang(ρ, k₋, k₊, γ) = γ * ρ * (1-ρ) * k₋ / (k₋ + k₊ + γ * ρ * fp(k₋,k₊)) #* Φ1(k₋, k₊)

    # oben boundary functions wang
    Wang_HD(β, k₋, k₊, γ) = γ * ρ_HD(β, γ) * (1-ρ_HD(β, γ)) * k₋ / (k₋ + k₊ + γ * ρ_HD(β, γ) * fp(k₋,k₊))
    Wang_HD_new(β, k₋, k₊, ϵ ) = ϵ * ρ_HD_new(β, ϵ, k₋, k₊) * (1-ρ_HD_new(β, ϵ, k₋, k₊)) * k₋ / (k₋ + k₊ + ϵ * ρ_HD_new(β, ϵ, k₋, k₊) * fp(k₋,k₊))
    
    Wang_LD(α, k₋, k₊, γ) = γ *   ρ_LD(α, k₋, k₊, γ)* (1-  ρ_LD(α, k₋, k₊, γ)) * k₋ / (k₋ + k₊ + γ *   ρ_LD(α, k₋, k₊, γ) * fp(k₋,k₊))
    
    ρ_HD(β, γ) = (1 - β/γ)
    ρ_LD(α, k₋, k₊, ϵ) = α * (k₊+k₋) / (ϵ * (k₋ - α*fp(k₋, k₊)))
    ρ_HD_new(β, ϵ, k₋, k₊) = (ϵ * k₋  - fa(k₋, k₊) * β* (k₋ + k₊)) / (β * ϵ *fp(k₋, k₊)*fa(k₋, k₊)+ ϵ * k₋)

    α_crit_fct(k₋, k₊, γ) = γ * ρ_max_fct(k₋,k₊, γ)/(1+ k₊/k₋)
    β_crit_fct(k₋, k₊, γ) = (1 - ρ_max_fct(k₋,k₊,γ)) * γ
    β_crit_fct_hat(k₋, k₊, γ, ρ) = β_crit_fct(k₋, k₊, γ) * fa_hat(k₋, k₊, ρ,γ)/ fa(k₋,k₊)
    k₊_crit_fct(k₋, α, γ) = (-2α + γ) / ((α/k₋)*(α/k₋+2)) #this is only for LD
    # Wang_k₊(k₊) = γ*ρ_max_fct(k₋,k₊,γ)*(1-ρ_max_fct(k₋,k₊,γ)) / ((1+f*τ)*(1+γf(k₊)*τ))
    
    #ρ_max = -(1+f*τ)/(γ*Φ(f, τ)*τ) + sqrt( ((1+f*τ) / (γ*Φ(f, τ)*τ))^2 + (1+f*τ)/(γ*Φ(f, τ)*τ))

    
    # α_crit = γ * ρ_max_fct(f, τ)/(1+ f*τ)
    active_fct(ρ, k₋, k₊) = ρ / (1+k₊/k₋)
    paused_fct(ρ, k₋, k₊) = ρ *(k₊/k₋)/(1+(k₊/k₋))
   
    # extended particle functions
    ρ_LD_ext(α, l_ribosome,ϵ) = α / (ϵ + α* (l_ribosome-1))
    J_LD_ext(α, ϵ, l_ribosome) = (α*(ϵ- α)) / (ϵ + α*(l_ribosome-1))
    α_crit_TASEP(l_ribosome,ϵ) = ϵ / ( 1 + sqrt(l_ribosome))

    ρ_HD_ext(β, l_ribosome,ϵ) = (1 - β/ϵ) / l_ribosome
    J_HD_ext(β, ϵ, l_ribosome) = (β*(1-β/ϵ)) / (1 - β/ϵ*(l_ribosome-1))
    β_crit_TASEP(l_ribosome,ϵ) = ϵ / ( 1 + sqrt(l_ribosome))
   
    J_ext_max(ϵ, l_ribosome) = ϵ  / (1+sqrt(l_ribosome))^2
    ρ_max_ext(l_ribosome) = 1 / (l_ribosome + sqrt(l_ribosome))

    J_ext(ϵ, ρ, l_ribosome) = (ϵ*ρ*(1-ρ*l_ribosome)) / (1 - ρ*(l_ribosome-1))
    J_ext_k₊(γ, l_ribosome, k₋, k₊, ρ) = k₋*γ*ρ*(1-ρ*l_ribosome) / ( (1-ρ*(l_ribosome-1)) * (k₋+ k₊ +γ*ρ*Φ(k₋,k₊)) ) 

    # Woosh functions
    pw_fct(ρ, L, k₋, k₊) = fa(k₋, k₊)^(L*ρ)

    P(d) = (1-Φ1(k₋, k₊))^d * Φ(k₋, k₊)
    J_whoosh(γ, ρ, L, k₋, k₊) = γ * ρ * (1-ρ) * P0(k₋, k₊, ρ,L)  + Wang(ρ, k₋, k₊, γ) * (1-P0(k₋, k₊, ρ,L) )
    Woosh_fraction(γ, ρ, L, k₋, k₊) = γ * ρ * (1-ρ) * pw_fct(ρ, L, k₋, k₊) / J_whoosh(γ, ρ, L, k₋, k₊) 
    Cluster_Wang_fraction(γ, ρ, L, k₋, k₊) = Wang(ρ, k₋, k₊, γ) * (1-pw_fct(ρ, L, k₋, k₊)) / J_whoosh(γ, ρ, L, k₋, k₊) 

    J_standard(ϵ, ρ) = ϵ*ρ*(1-ρ)
    ρ_LD_standard(α, ϵ) = α / ϵ
    ρ_LD_standard_finitesize( k₋, k₊, L, ρ) = ρ* P0(k₋, k₊, ρ,L)
    

    J_LD_standard_finitesize( k₋, k₊, L,ϵ, ρ) = P0(k₋, k₊, ρ,L) * J_standard(ϵ, ρ)
    J_MC_standard_finitesize( k₋, k₊, L,ϵ) = P0(k₋, k₊, 0.5,L) * J_standard(ϵ, 0.5)
    
    J_pTASEP_LD_standard_finitesize(α, k₋, k₊, ϵ)= (1-P0(k₋, k₊, α/ϵ, L)) * Wang_LD(α, k₋, k₊, ϵ)
    J_pTASEP_MC_standard_finitesize(k₋, k₊, ρ, ϵ)= (1-P0(k₋, k₊, 0.5, L)) * Wang(ρ, k₋, k₊, ϵ)

    ψ(ρ,k₋, k₊) = ρ *(k₊)/(k₋+k₊)
    ν(ρ, k₋, k₊, γ ) = ρ - Wang(ρ, k₋, k₊, γ) / γ - ψ(ρ,k₋, k₊)
    μ(ρ, k₋, k₊, γ) = Wang(ρ, k₋, k₊, γ) / γ
    
    greulich(ρ, k₋, k₊, γ) = γ*ρ*(1-ρ) / (1 + (γ/k₋) * ρ * (k₊ / (k₊ + k₋) ) )

end


 ############################################################### Lorenzo functions ######################################################################
function av_current_half(ρ,L,γ,t)
    if t<ρ*L
        return γ/L*(t/6)
    end
    
    if t>=ρ*L
        if t<L/(4*ρ)
            return (γ*ρ^2*L/6+(γ*ρ/3)*(ρ*L+3*t-4*sqrt(ρ*L*t)))/t
        end
        if t>=L/(4*ρ)
            return (ρ^2*L/6+(ρ/3)*(ρ*L+3*L/(4*ρ)-2*L)+ρ*(1-ρ)*t-L/4*(1-ρ)-(L^2/48)*(4*ρ/L-1/t))*(γ/t)
        end
    end

end

function find_crit_value_num(α, k₋, k₊_list, γ)d
    for k₊ in k₊_list

        J_W_LD = Wang_LD(α, k₋, k₊, γ)
        J_W_MC = Wang(ρ_max_fct(k₋,k₊, γ), k₋, k₊, γ)
        if  0.99 <= J_W_LD/ J_W_MC <= 1.01
            return k₊
        end
    end
end

function obc_current_wang(ρ, k₊, k₋, γ, α,kp_crit)
    if k₊ < kp_crit
        return Wang_LD(α, k₋, k₊, γ)
    else
        return Wang.(ρ, k₋, k₊, γ)
    end
end


function trans_cur(ρ,L,γ,t)
    if ρ<=0.5
        return istant_current_half(ρ,L,γ,t)
    end
    return istant_current_half(1-ρ,L,γ,t)
end

function av_trans_cur(ρ,L,γ,t)
    if ρ<=0.5
        return av_current_half(ρ,L,γ,t)
    end
    return av_current_half(1-ρ,L,γ,t)
end


function J_transient_w(ρ,kp,km,L)
    #new_prob=(ρ*L-1)*kp/((ρ*L-1)*kp+km)*(1-fa)*fa^(ρL-1)
    fa = 1- kp/(km+kp)
    D = 1/(1-fa)
    P₀=fa^(ρ*L)
    if ρ==0
        return 0
    else
    return P₀*av_trans_cur(ρ,L,1,1/(ρ*L*kp))         
    end
end

function J_transient_c(ρ,kp,km,L)
    #new_prob=(ρ*L-1)*kp/((ρ*L-1)*kp+km)*(1-fa)*fa^(ρL-1)
    fa = 1- kp/(km+kp)
    D = 1/(1-fa)
    P₀=fa^(ρ*L)
    if ρ==0
        return 0
    else
    return (1-P₀)*km*(1-ρ)*ρ*L*D/(ρ*L+D-1)                #-(ρ*L-1)*kp/((ρ*L-1)*kp+km)*(1-fa)*fa^(ρ*L-1)
    end
end


function J_transient(ρ,kp,km,L)
    #new_prob=(ρ*L-1)*kp/((ρ*L-1)*kp+km)*(1-fa)*fa^(ρL-1)
    fa = 1- kp/(km+kp)
    D = km /kp * (1 - (km / (km + kp))^(ρ * L))
    P₀=fa^(ρ*L)
    if ρ==0
        return 0
    else
    return P₀*av_trans_cur(ρ,L,1,1/(ρ*L*kp))+ (1-P₀)*km*(1-ρ)*ρ*L*D/(ρ*L+D-1)                #-(ρ*L-1)*kp/((ρ*L-1)*kp+km)*(1-fa)*fa^(ρ*L-1)
    end
end




############################################################################ weights ######################################################################


function distance_paused_particles(k₋, k₊, N)
    d = (k₋ / k₊) * (1 - (k₋ / (k₋ + k₊))^N)
    return d
end

function d(N, kp, km)
    r = km / (km + kp)
    return (km / kp) * (1 - r^N)
end

function time_clogged_OBC_non_iter(k₋, k₊, l_ribosome,avg_pausing_pos)
    N = avg_pausing_pos / l_ribosome
    d = distance_paused_particles(k₋, k₊, N)
    #d = k₋ / k₊ + 1
    τ_c = N/(d+1) * 1/k₋

    return τ_c

end

function time_clogged_OBC_non_iter_inf(k₋, k₊, l_ribosome,avg_pausing_pos)
    N = avg_pausing_pos / l_ribosome
    #d = distance_paused_particles(k₋, k₊, N)
    d = k₋ / k₊ + 1
    τ_c = N/(d+1) * 1/k₋

    return τ_c

end

function weight(L,kp,km,l_ribosome, init, elong,avg_pausing_pos)
    
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        ρ = ρ_max_ext(l_ribosome)
    else
        ρ = ρ_LD_ext(init, l_ribosome,elong)
    end
    w1 = analyitcal_total_expected_pausing_time_extended(init,elong,L,kp,l_ribosome)
    w2 = time_clogged_OBC_non_iter(km, kp, l_ribosome, avg_pausing_pos)
    #w2 = cluster_lifespan_light_extended(L,kp,km,l_ribosome)
    return w1 / (w2 + w1)
    #fp=kp/(kp+km)
    #if kp != 0
    #    w1 = analyitcal_total_expected_pausing_time_extended(init,elong,L,kp,l_ribosome)
    #    w2 = time_clogged_OBC_non_iter(km, kp, l_ribosome, avg_pausing_pos)
        #w2 = cluster_lifespan_light_extended(L,kp,km,l_ribosome)
    #    return w1 / (w2 + w1)
    #else
    #    return 1
    #end
end


function time_in_paused(init, l_ribosome, elong , L, kp, km)
    
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        ρ0 = ρ_max_ext(l_ribosome)
    else
        ρ0 = ρ_LD_ext(init, l_ribosome,elong)
    end

    avg_pausing_pos = analyitcal_total_expected_pausing_position_extended(init,elong,L,kp,l_ribosome)
    w2 = time_clogged_OBC_non_iter(km, kp, l_ribosome, avg_pausing_pos)

    if kp == 0
        return 0
    else
        return w2
    end

end

function weight_inf(L,kp,km,l_ribosome, init, elong,avg_pausing_pos)
    
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        ρ = ρ_max_ext(l_ribosome)
    else
        ρ = ρ_LD_ext(init, l_ribosome,elong)
    end
    #fp=kp/(kp+km)
    if kp != 0
        w1 = analyitcal_total_expected_pausing_time_extended(init,elong,L,kp,l_ribosome)
        w2 = time_clogged_OBC_non_iter_inf(km, kp, l_ribosome, avg_pausing_pos)
        #w2 = cluster_lifespan_light_extended(L,kp,km,l_ribosome)
        return w1 / (w2 + w1)
    else
        return 1
    end
end

function cluster_dissolution_time(L, kp, km, l_ribosome, avg_pausing_pos)
    # If there is no pausing, we return 0.
    if kp == 0
        return 0.0
    end

    total_dissolution_time = 0.0

    # Determine the number of simulation runs based on L.
    number_of_runs = L <= 500 ? 500 : 100

    # Compute maximum number of cluster elements to consider.
    max_nc = round(Int, 2 * avg_pausing_pos / l_ribosome)

    
    # Loop over simulation runs.
    for run in 1:number_of_runs
        run_time = 0.0

        # For each cluster element index (interpreted as the threshold that must be met)
        for cluster in 1:max_nc
            d = 0.0           # current accumulated distance
            c_time = 0.0      # cumulative time for this cluster element

            # Continue accumulating time until the distance d reaches or exceeds the threshold 'cluster'
            while d < cluster
                t_added = rand(Exponential(1/km))
                c_time += t_added

                # p depends on the current cumulative time
                p = kp/(km+kp) * (1 - exp(-(kp+km)*c_time))
                
                # Increase d by 1 plus the mean of a NegativeBinomial(1, p)
                # (For NegativeBinomial(1,p), mean = (1-p)/p if p>0)
                d += mean(NegativeBinomial(1, p)) + 1
            end

            # Add the dissolution time for this cluster element
            run_time += c_time
        end

        # Average over the cluster elements for this run
        total_dissolution_time += run_time / max_nc
    end

    # Return the average dissolution time over all simulation runs.
    return total_dissolution_time / number_of_runs
end

function weight_itr(L,kp,km,l_ribosome, init, elong, avg_pausing_pos)
    
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        ρ = ρ_max_ext(l_ribosome)
    else
        ρ = ρ_LD_ext(init, l_ribosome,elong)
    end

    #fp=kp/(kp+km)
    if kp != 0
        w1=analyitcal_total_expected_pausing_time_extended(init,elong,L,kp,l_ribosome)
        w2=cluster_dissolution_time(L, kp, km, l_ribosome, avg_pausing_pos)
        return w1 / (w2 + w1)
    else
        return 1
    end
end

################################################################################ density functions ##########################################################

function closed_form_density(N_max, kp, km, L, l_ribosome) 
    #analyitcal_function for the density in the paused state
    B = km/kp + 1.0
    if B >= L/l_ribosome
        B = L/l_ribosome
    end
    
    density = ((2 * N_max + 1)/6 + B/2) / L
    return density
end

function total_density_inf(init, l_ribosome, elong , L, kp, km)
    
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        ρ0 = ρ_max_ext(l_ribosome)
    else
        ρ0 = ρ_LD_ext(init, l_ribosome,elong)
    end

    avg_pausing_pos = analyitcal_total_expected_pausing_position_extended(init,elong,L,kp,l_ribosome)
    w1 = weight(L,kp,km,l_ribosome, init, elong, avg_pausing_pos)
    N_max = round(2*avg_pausing_pos / l_ribosome)
    if kp == 0
        ρC = 0
    else
        ρC = closed_form_density(N_max, kp, km, L,l_ribosome)
    end

    return w1 * ρ0 + (1-w1) * ρC
end

function unpaused_density_inf(init, l_ribosome, elong , L, kp, km)
    
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        ρ0 = ρ_max_ext(l_ribosome)
    else
        ρ0 = ρ_LD_ext(init, l_ribosome,elong)
    end

    avg_pausing_pos = analyitcal_total_expected_pausing_position_extended(init,elong,L,kp,l_ribosome)
    w1 = weight(L,kp,km,l_ribosome, init, elong, avg_pausing_pos)
    N_max = round(2*avg_pausing_pos / l_ribosome)
    if kp == 0
        ρC = 0
    else
        ρC = closed_form_density(N_max, kp, km, L,l_ribosome)
    end

    return w1 * ρ0
end

function paused_density_inf(init, l_ribosome, elong , L, kp, km)
    
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        ρ0 = ρ_max_ext(l_ribosome)
    else
        ρ0 = ρ_LD_ext(init, l_ribosome,elong)
    end

    avg_pausing_pos = analyitcal_total_expected_pausing_position_extended(init,elong,L,kp,l_ribosome)
    w1 = weight(L,kp,km,l_ribosome, init, elong, avg_pausing_pos)
    N_max = round(2*avg_pausing_pos / l_ribosome)
    if kp == 0
        ρC = 0
    else
        ρC = closed_form_density(N_max, kp, km, L,l_ribosome)
    end

    return (1-w1) * ρC
end

function paused_density_inf_not_weighted(init, l_ribosome, elong , L, kp, km)
    
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        ρ0 = ρ_max_ext(l_ribosome)
    else
        ρ0 = ρ_LD_ext(init, l_ribosome,elong)
    end

    avg_pausing_pos = analyitcal_total_expected_pausing_position_extended(init,elong,L,kp,l_ribosome)
    #w1 = weight_inf(L,kp,km,l_ribosome, init, elong, avg_pausing_pos)
    w1 = 1 - weight(L,kp,km,l_ribosome, init, elong, avg_pausing_pos)
    N_max = round(2*avg_pausing_pos / l_ribosome)
    #if kp == 0
    #    ρC = 0
    #else
        ρC = closed_form_density(N_max, kp, km, L,l_ribosome)
    #end
    #return w1 *ρC
    return ρC
end

function paused_density_inf_weighted(init, l_ribosome, elong , L, kp, km)
    
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        ρ0 = ρ_max_ext(l_ribosome)
    else
        ρ0 = ρ_LD_ext(init, l_ribosome,elong)
    end

    avg_pausing_pos = analyitcal_total_expected_pausing_position_extended(init,elong,L,kp,l_ribosome)
    #w1 = weight_inf(L,kp,km,l_ribosome, init, elong, avg_pausing_pos)
    w1 = 1 - weight(L,kp,km,l_ribosome, init, elong, avg_pausing_pos)
    N_max = round(2*avg_pausing_pos / l_ribosome)
    #if kp == 0
    #    ρC = 0
    #else
        ρC = closed_form_density(N_max, kp, km, L,l_ribosome)
    #end
    return w1 * ρC
    #return ρC
end

function cluster_average_density(L, kp, km, l_ribosome, avg_pausing_pos)
    # If there is no pausing, we return 0.
    #if kp == 0
    #    return 0.0
    #end

    # Determine the number of simulation runs based on L.
    number_of_runs = L <= 200 ? 1000 : 100

    # Compute the maximum number of cluster elements to consider.
    max_nc = round(Int, 2 * avg_pausing_pos / l_ribosome)

    total_density_runs = 0.0  # Accumulates the average density from each run

    for run in 1:number_of_runs
        run_total_density = 0.0   # Accumulates density (time-integrated) over all clusters in this run
        run_total_time_steps = 0  # Accumulates the total number of time steps over all clusters

        # Loop over each cluster (initial cluster sizes from 1 to max_nc)
        for cluster_index in 1:max_nc
            N_c = cluster_index         # Initial number of particles in the cluster
            cluster_density_sum = 0.0   # Accumulates density over time for this cluster
            cluster_time_steps = 0      # Counts time steps for this cluster simulation
            c_time = 0.0                # Cumulative time for this cluster simulation

            # Simulate the dissolution until all particles in the cluster have unpaused.
            while N_c > 0
                # Calculate the density at the current time step.
                current_density = N_c / L
                cluster_density_sum += current_density
                cluster_time_steps += 1

                # Increment time by drawing a random waiting time.
                t_added = rand(Exponential(1/km))
                c_time += t_added

                # Calculate the probability of an unpausing event at the current time.
                p = kp / (km + kp) * (1 - exp(- (kp + km) * c_time))
                
                # Determine the expected number of unpausing events.
                # For a NegativeBinomial(1, p), mean = (1-p)/p, so we add 1.
                d = mean(NegativeBinomial(1, p))
                expected_unpauses = d + 1
                
                #expected_unpauses = distance_paused_particles_new(km, kp, cluster_index) + 1

                # Update the number of particles remaining in the cluster.
                N_c -= expected_unpauses
                if N_c < 0
                    N_c = 0
                end
            end

            # Scale and accumulate the cluster's density.
            # (Multiplying by 1/km corresponds to scaling by the mean waiting time.)
            run_total_density += cluster_density_sum * (1/km)
            run_total_time_steps += cluster_time_steps
        end

        # Calculate the average density for this run.
        # The factor 1/km cancels if applied both in accumulation and in averaging.
        average_density_run = run_total_density / (run_total_time_steps * (1/km))
        total_density_runs += average_density_run
    end

    # Return the average density over all simulation runs.
    return total_density_runs / number_of_runs
end

function ρ0_ext_itr(init, l_ribosome,elong,L,kp,km)

    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        ρ = ρ_max_ext(l_ribosome)
    else
        ρ = ρ_LD_ext(init, l_ribosome,elong)
    end
    #densC = dens_c(km, kp, L, l_ribosome) / L 
    avg_pausing_pos = analyitcal_total_expected_pausing_position_extended(init,elong,L,kp,l_ribosome)
    w1 = weight_itr(L,kp,km,l_ribosome, init, elong, avg_pausing_pos)
    #w1 = weight_itr(L,kp,km,l_ribosome, dens0)

    return w1 * ρ 

end

function ρc_ext_itr(init, l_ribosome,elong,L,kp,km)

    avg_pausing_pos = analyitcal_total_expected_pausing_position_extended(init,elong,L,kp,l_ribosome)
    N_max = round(2*avg_pausing_pos / l_ribosome)
    densC = numeric_paused_density(N_max, L, km, kp, B_of, l_ribosome)#cluster_average_density(L, kp, km, l_ribosome, avg_pausing_pos)
    #densC = avg_pausing_pos/(2*l_ribosome*L)
    w2 = 1 - weight(L,kp,km,l_ribosome, init,elong ,avg_pausing_pos)

    return w2  * densC 

end

function ρc_ext_itr_not_weighted(init, l_ribosome,elong,L,kp,km)

    avg_pausing_pos = analyitcal_total_expected_pausing_position_extended(init,elong,L,kp,l_ribosome)
    N_max = round(2*avg_pausing_pos / l_ribosome)
    #densC = avg_pausing_pos/(2*l_ribosome*L)
    #w2 = 1 - weight_itr(L,kp,km,l_ribosome, init,elong ,avg_pausing_pos)
    #w2 = 1 - weight(L,kp,km,l_ribosome, init,elong ,avg_pausing_pos)
    #if kp == 0
    #    densC = 0
    #else
        densC = numeric_paused_density(N_max, L, km, kp, B_of, l_ribosome)#cluster_average_density(L, kp, km, l_ribosome, avg_pausing_pos)

    #end
    return densC 

end

function ρ_ext_OBC_itr(α, l_ribosome,ϵ,L,kp,km)
    return ρ0_ext_itr(α, l_ribosome,ϵ,L,kp,km) + ρc_ext_itr(α, l_ribosome,ϵ,L,kp,km)
end

"""
    numeric_paused_density(Cmax, L, B_of)

"""
function numeric_paused_density(N_max, L, km, kp, B_of, l_ribosome)
    total_weight     = 0.0
    weighted_density = 0.0
    for Ni in 1:N_max
        # block‐size for this Ni
        Bi = B_of(km, kp, Ni)

        if Bi >= L/l_ribosome
            Bi = L/l_ribosome
        end
        # number of blocks
        nB = ceil(Int, Ni / Bi)

        # mean occupied sites per block:
        #   <N_c> = (nB*Ni - Bi*(nB-1)*nB/2) / nB
        N_c = (Ni + Bi) /2

        total_weight     += nB
        weighted_density += nB * (N_c / L)
    end

    return weighted_density / total_weight
end

B_of(km,kp,Ni) = distance_paused_particles(km, kp, Ni) + 1

function total_density_numerics(init, l_ribosome, elong , L, kp, km)
    
    if init/elong >= 1/ (1 + sqrt(l_ribosome))
        ρ0 = ρ_max_ext(l_ribosome)
    else
        ρ0 = ρ_LD_ext(init, l_ribosome,elong)
    end

    B_of(km,kp,Ni) = distance_paused_particles(km, kp, Ni) + 1
    
    avg_pausing_pos = analyitcal_total_expected_pausing_position_extended(init,elong,L,kp,l_ribosome)
    w1 = weight(L,kp,km,l_ribosome, init, elong, avg_pausing_pos)
    
    N_max = round(2*avg_pausing_pos / l_ribosome)
    
    if kp == 0
        ρC = 0
    else
        ρC = numeric_paused_density(N_max, L, km, kp, B_of, l_ribosome)
    end

    return w1 * ρ0 + (1-w1) * ρC
end


#################################################################################### current ##########################################################

function Jc_extended(L,kp,km,l_ribosome, init, elong)
    if kp == 0.0
        return 0
    else
        avg_pausing_pos = analyitcal_total_expected_pausing_position(init,elong,L,kp)
        τ1 = time_clogged_OBC_non_iter(km, kp, l_ribosome, avg_pausing_pos)
        #τ1 = cluster_lifespan_light_extended(L,kp,km,l_ribosome)
        w = weight(L,kp,km,l_ribosome, init, elong, avg_pausing_pos)
        J = (1 - w) * avg_pausing_pos/(τ1*l_ribosome)
        #J =  L/(2*τ1*l_ribosome)
        return J
    end
end

function J0_extended(L,kp,km,l_ribosome, init, elong)
    avg_pausing_pos = analyitcal_total_expected_pausing_position(init,elong,L,kp)
    w = weight(L,kp,km,l_ribosome, init,elong,avg_pausing_pos)
    J = w * J_LD_ext(init, elong, l_ribosome)
    return J
end


function current_OBC(L,kp,km,l_ribosome,init, elong)
    return J0_extended(L,kp,km,l_ribosome, init, elong) + Jc_extended(L,kp,km,l_ribosome, init, elong)
end



function Jc_extended_itr(L,kp,km,l_ribosome, init, elong)
    if kp == 0.0
        return 0
    else
        avg_pausing_pos = analyitcal_total_expected_pausing_position(init,elong,L,kp)
        ρ = ρ_LD_ext(init, l_ribosome, elong)
        #τ1 = time_clogged_OBC_non_iter(km, kp, L, l_ribosome)
        τ1 = cluster_dissolution_time(L, kp, km, l_ribosome, avg_pausing_pos)
        w = weight_itr(L,kp,km,l_ribosome, init, elong, avg_pausing_pos)
    
        J = (1 - w) * L/(2*τ1*l_ribosome)
        #J =  L/(2*τ1*l_ribosomeƒ)
        return J
    end
end


function J0_extended_itr(L,kp,km,l_ribosome, init, elong)
    avg_pausing_pos = analyitcal_total_expected_pausing_position(init,elong,L,kp)
    w = weight_itr(L,kp,km,l_ribosome, init,elong,avg_pausing_pos)
    J = w * J_LD_ext(init, elong, l_ribosome)
    return J
end


function current_OBC_itr(L,kp,km,l_ribosome, α, ϵ)
    return  J0_extended_itr(L,kp,km,l_ribosome, α, ϵ) + Jc_extended_itr(L,kp,km,l_ribosome, α, ϵ)
end


###################################################################### instantaneousDensityProfile ##########################################################



function instantaneousDensityProfile_ext(L::Int,l_ribosome::Int, alpha::Float64, epsilon::Float64, t::Float64) 
    dens = Vector{Float64}(undef, L)
    # Shock speed from left boundary
    vshock = epsilon - alpha  # rankine-hugoniot for rho=alpha/epsilon -> 0

    # Time for the shock to reach the right boundary
    t_star = L / vshock  

    # Boundary density behind the shock
    rho_behind = alpha /(elong+alpha*(l_ribosome-1))

    if t < 0
        # negative time doesn't make sense here, so return zeros
        fill!(dens, 0.0)
        return dens
    end

    if t < t_star
        # The shock is at position x_s(t) = vshock * t
        x_shock = vshock * t

        # Fill sites i < x_shock with alpha/epsilon, else 0
        for i in 1:L
            if i < x_shock
                dens[i] = rho_behind
            else
                dens[i] = 0.0
            end
        end
    else
        # Shock has already traversed the entire lattice
        # => the whole domain is behind the shock with density alpha/epsilon
        fill!(dens, rho_behind)
    end

    return dens
end


"""
    totalParticleNumber(L, alpha, epsilon, t)

Compute N(t), the total number of particles on the lattice at time t,
using the same short-time shock-based TASEP approximation.

Arguments:
- L::Int           : Number of lattice sites (length)
- alpha::Float64   : Injection rate [1/time]
- epsilon::Float64 : Hopping rate [1/time]
- t::Float64       : Time

Returns:
- A Float64 representing N(t).
  
Short-time formula:
- For t < t_star = L/(ε - α):
    N(t) = (α/ε)*(ε - α)* t
  i.e. α(1-α/ε)·t
- For t ≥ t_star:
    N(t) = (α/ε)*L

Same assumptions as above.
"""
function totalParticleNumber(L::Int, alpha::Float64, epsilon::Float64, t::Float64)
    # shock speed
    vshock = epsilon - alpha
    t_star = L / vshock
    rho_behind = alpha / epsilon

    if t < 0
        return 0.0
    elseif t < t_star
        # N(t) = integral( rho behind ) from x=0.. x=(vshock t)
        #      = (vshock t)*(alpha/epsilon) = (epsilon - alpha)*t * alpha/epsilon
        return (epsilon - alpha)*t * (alpha / epsilon) 
    else
        # entire domain is behind the shock => N(t) = alpha/epsilon * L
        return rho_behind * L
    end
end


##################### survival probability distribution ################################
# Calculate effective arrival rate R and characteristic time t_star
function f_exponential(t, init, elong, kp, L)
    t_star = L / (elong - init)
    lambda_out = (init/ elong) * L * kp
    return lambda_out*exp(-lambda_out*(t-t_star*0.5))
end

function f_rayleigh(t, init, elong, kp)
    R = init * (1 - (init /elong))
    return R * kp * t * exp(-0.5 * R * kp * t^2)
end

# Combined PDF
function f_combined(t, init, elong, kp, L)
    t_star = L / (elong - init)
    if t <= t_star
        return f_rayleigh(t, init, elong, kp)
    else
        return f_exponential(t, init, elong, kp, L)
    end
end


##################### Expected time from unpaused state ################################

function analytical_expected_value_exp_time(init, elong, L, kp)
    #this functions calculates the analyrical expected time until a pausing event occurs for ell = 1
    t_star = L / (elong - init)
    lambda = (init/ elong) * L * kp
    return exp(-0.5*lambda*t_star) *(t_star + 1/ lambda)
end

function analytical_expected_value_raiyleigh_time(init, elong, L, kp)

    #this functions calculates the analyrical expected time until a pausing event occurs for ell = 1 before domain wall exited the lattice
    t_star = L / (elong - init)
    R = init * (1 - (init /elong))
    a = 0.5*R*kp
    return sqrt(pi)/(2*sqrt(a))*erf(sqrt(a)*t_star) - t_star*exp(-a*t_star^2)
end


function analytical_expected_value_exp_time_extended(init,elong,L,kp,l_ribosome)
    if init/elong >= 1/(1+sqrt(l_ribosome))
        t_star = L *(1+sqrt(l_ribosome))/(elong*l_ribosome)
        lambda = L/(l_ribosome+sqrt(l_ribosome)) * kp
    else
        t_star = L / (elong - init)
        lambda = (init /(elong + init*(l_ribosome-1)))*L *kp
    end
    return exp(-0.5*lambda*t_star) *(t_star + 1/ lambda)
end

function analytical_expected_value_exp_time_extended2(init,elong,L,kp,l_ribosome)
    if init/elong >= 1/(1+sqrt(l_ribosome))
        t_star = L *(1+sqrt(l_ribosome))/(elong*l_ribosome)
        lambda = L/(l_ribosome+sqrt(l_ribosome)) * kp
    else
        t_star = L / (elong - init)
        lambda = (init /(elong + init*(l_ribosome-1)))*L *kp
        R = init * (elong - init ) / (elong + init*(l_ribosome - 1))
        a = 0.5*R*kp
    end
    return exp(-a*t_star^2) *(t_star + 1/ lambda)#exp(-lambda*t_star) *(t_star + 1/ lambda) #
end

function analytical_expected_value_raiyleigh_time_extended(init, elong, L, kp, l_ribosome)
    if init/elong >= 1/(1+sqrt(l_ribosome))
        t_star = L *(1+sqrt(l_ribosome))/(elong*l_ribosome)
        R = elong /(1 + sqrt(l_ribosome))^2
        a = 0.5*R*kp
    else
        t_star = L / (elong - init)
        R = init * (elong - init ) / (elong + init*(l_ribosome - 1))
        a = 0.5*R*kp
    end

    return sqrt(pi)/(2*sqrt(a)) * erf(sqrt(a)*t_star) - t_star*exp(-a*t_star^2)
end

function analyitcal_total_expected_pausing_time(init,elong,L,kp)
    E1 = analytical_expected_value_exp_time(init, elong, L, kp)
    E2 = analytical_expected_value_raiyleigh_time(init, elong, L, kp)
    return E1 + E2
end

function analyitcal_total_expected_pausing_time_extended(init,elong,L,kp,l_ribosome)
    E1 = analytical_expected_value_exp_time_extended(init,elong,L,kp,l_ribosome)
    E2 = analytical_expected_value_raiyleigh_time_extended(init, elong, L, kp, l_ribosome)
    return E1 + E2
end


##################### Expected time from unpaused state ################################

function analytical_expected_value_exp_position(init, elong, L, kp)
    t_star = L / (elong - init)
    lambda = (init/ elong) * L * kp
    return exp(-0.5*lambda*t_star) * L / 2
end

function analytical_expected_value_raiyleigh_position(init, elong, L, kp)
    return 0.5*(elong - init) * analytical_expected_value_raiyleigh_time(init, elong, L, kp)
end


function analyitcal_total_expected_pausing_position(init,elong,L,kp)
    E1 = analytical_expected_value_exp_position(init, elong, L, kp)
    E2 = analytical_expected_value_raiyleigh_position(init, elong, L, kp)
    return E1 + E2
end


function analytical_expected_value_exp_position_ext(init, elong, L, kp, l_ribosome)
    if init/elong >= 1/(1+sqrt(l_ribosome))
        t_star = L *(1+sqrt(l_ribosome))/(elong*l_ribosome)
        lambda = L/(l_ribosome+sqrt(l_ribosome)) * kp
    else
        t_star = L / (elong - init)
        lambda = (init /(elong + init*(l_ribosome-1)))*L *kp
    end
    return exp(-0.5*lambda*t_star) * L / 2
end

function analytical_expected_value_raiyleigh_position_extended(init, elong, L, kp,l_ribosome)
    if init/elong >= 1/(1+sqrt(l_ribosome))
        s = elong*sqrt(l_ribosome) /(1 + sqrt(l_ribosome))
    else
        s = 0.5*(elong - init)
    end
    return  s* analytical_expected_value_raiyleigh_time_extended(init, elong, L, kp,l_ribosome)
end

function analyitcal_total_expected_pausing_position_extended(init,elong,L,kp,l_ribosome)
    E1 = analytical_expected_value_raiyleigh_position_extended(init, elong, L, kp,l_ribosome)
    E2 = analytical_expected_value_exp_position_ext(init, elong, L, kp, l_ribosome)
    return E1 + E2
end

function pdf_pausing_position_extended(x, init,elong, kp, L, l_ribosome)
    v = elong - init
    a = J_LD_ext(init, elong, l_ribosome) * kp / 2
    if x < L
        return (2 * a * x / v^2) * exp(-a * x^2 / v^2)
    else
        return exp(-a * L^2 / v^2)  # P_after at x = L
    end
end


