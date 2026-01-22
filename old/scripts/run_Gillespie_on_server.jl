using Random
using Distributions
using Distributed
using JLD
#using Colors
#using BenchmarkTools
using Pkg

number_of_workers = 2
addprocs(number_of_workers)

@everywhere include("Gillespie_obc.jl")
@everywhere include("Basic_functions_obc.jl")
@everywhere include("analyitcal_functions.jl")


pmap(_ -> precompile_task(), 1:nworkers())

function save_data_to_dict(outputs, data_dict, simulation_results, listby)
    for (index, result) in enumerate(simulation_results)
        c = listby[index]
        for (j, output) in enumerate(outputs)
            data_dict["Results"][output][c] = result[j]
        end
    end
    return data_dict
end

L, l_ribosome, track_site = 250, 1, 1
init, term, elong = 1, 1, 1
deg_t = deg_g = 0
k₋ = 0.0001
k₊ = 0.00001 

#α_crit = α_crit_hat_fct(k₋, k₊, elong, ρ_max_fct(k₋, k₊, 1))
α_list = 10 .^ collect(-6:0.25:0)
#β_crit = β_crit_fct(k₋, k₊, 1) * fa_hat(k₋, k₊, ρ_max_fct(k₋, k₊, 1),1)/ fa(k₋,k₊)
β_list = collect(range(start = 0.01, stop = 0.6, length = 36))


x = 1
run_time, starting_t, delta_t = 4 * 10^x, 1*10^x, 3*10^x


function make_phase_diagram(L, l_ribosome, track_site, deg_t, deg_g, α_list, β_list, elong, k₋, k₊, run_time, starting_t, delta_t)
    
    phase_diagram_current = Matrix{Any}(undef,length(α_list), length(β_list))
    phase_diagram_density = Matrix{Any}(undef, length(α_list), length(β_list))
    phase_diagram_paused_distr = Matrix{Any}(undef, length(α_list), length(β_list))
    phase_diagram_number_distr = Matrix{Any}(undef, length(α_list), length(β_list))
    phase_diagram_length_distr = Matrix{Any}(undef, length(α_list), length(β_list))
    phase_diagram_time = Matrix{Any}(undef, length(α_list), length(β_list))
    phase_diagram_rates = Matrix{Any}(undef, length(α_list), length(β_list))

    i = 0
    for init in α_list

        #if init > α_list[Int64(length(α_list)/2)]
        #    run_time /= 2
        #    starting_t /= 2
        #    delta_t /= 2
        #end 

        results = pmap(
            β -> Gillespie_obc(L, l_ribosome, track_site, deg_t, deg_g, init, β, elong, k₋, k₊, run_time, starting_t, delta_t; kymo=false),
            β_list)

        i += 1    
        for (index, result) in enumerate(results)

            phase_diagram_current[i,index] = result[1]
            phase_diagram_density[i,index] = result[2]
            #phase_diagram_paused_distr[i,index] = result[3]
            #phase_diagram_number_distr[i,index] = result[4]
            #phase_diagram_length_distr[i,index] = result[6]
            phase_diagram_time[i,index] = result[3]
        end 

    end
    #return phase_diagram_current, phase_diagram_density, phase_diagram_paused_distr, phase_diagram_number_distr, phase_diagram_length_distr, phase_diagram_time
    return phase_diagram_current, phase_diagram_density, phase_diagram_time
end

@time PD_current, PD_density,  PD_time = make_phase_diagram(L, l_ribosome, track_site, deg_t, deg_g, α_list[1:2], β_list[1:2], elong, k₋, k₊, run_time, starting_t, delta_t)

x = 4
run_time, starting_t, delta_t = 12 *10^x, 4*10^x, 8*10^x 

@time PD_current, PD_density, PD_time = make_phase_diagram(L, l_ribosome, track_site, deg_t, deg_g, α_list, β_list, elong, k₋, k₊, run_time, starting_t, delta_t)




@save "PD_current_shifted.jld" PD_current
@save "PD_density_shifted.jld" PD_density
#@save "PD_paused_distr_shifted.jld" PD_paused_distr
#@save "PD_number_distr_shifted.jld" PD_number_distr
#@save "PD_length_distr_shifted.jld" PD_length_distr
@save "PD_time_shifted.jld" PD_time


#run_time, starting_t, delta_t = 40 * 10^x, 10*10^x, 30*10^x
#delta_t / (1/k₊)
#L_list = [500, 1000, 1500, 2000]
