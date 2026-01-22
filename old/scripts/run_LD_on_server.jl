using Random
using Plots
using Distributions
using Distributed
using JLD
using QuadGK         # for numerical integration
using SpecialFunctions  # for erf and sqrt(pi)
using StatsBase
using SavitzkyGolay
nworkers()
number_of_workers = Sys.CPU_THREADS + 2 
addprocs(number_of_workers)
rmprocs(1:number_of_workers)

@everywhere include("Gillespie_obc.jl")
@everywhere include("Basic_functions_obc.jl")
#include("density_profile_domain_wall.jl")
@everywhere include("analyitcal_functions.jl")
@everywhere include("degradation_functions.jl")
include("Gillespie_obc.jl")

function save_data_to_dict(outputs, data_dict, simulation_results, listby)
    for (index, result) in enumerate(simulation_results)
        c = listby[index]
        for (j, output) in enumerate(outputs)
            data_dict["Results"][output][c] = result[j]
        end
    end
    return data_dict
end


outputs = [
    "current", 
    "density", 
    "paused_part_distribution", 
    "cluster_length_distribution", 
    "cluster_number_distribution_unp",
    "cluster_number_distribution_p",
    "total_time",
    ]


L, l_ribosome, track_site = 850, 10, 4
ell = l_ribosome
ϵ = 20.0
init, term, elong = 0.2, ϵ, ϵ
deg_t = deg_g = 1/120
conc = range(0.1, stop = 12, length = 24)
km = 0.084 / 60
kp =  0.034 / 60 
kp_list = 10 .^ collect(-5:0.1:-1)
conc = kp_list ./ kp

kp_list = 10 .^ collect(-5:0.1:-1)
kp_list = collect(kp .* conc)
conc_ana = collect(-5:0.1:-1)

#KD = 0.5
#k₋ = KD * k₊
#RNAP pausing
#k₋ = 0.5#/ ϵ
#k₊ = 0.25#/ ϵ
#=
#tet
KD = 5 #μM the range given in papers is 0.5 - 20 μM
k₊ = 0.285 / ϵ # μM per second measured by tritton et.al 
k₋ = KD * k₊ / ϵ
=#

begin
    starting_t, delta_t = 10^2, 4 * 10^6
    run_time = starting_t + delta_t

    @time results850_tau120_paper_log = pmap(
    kp -> Gillespie_obc(L, l_ribosome, track_site, deg_t, deg_g, init, term, elong, km, kp, run_time, starting_t, delta_t; kymo=false),
    kp_list)

    saving_dict = Dict(
        "Input" => Dict(
            "L" => L,
            "l_ribosome" => l_ribosome,
            "track_site" => track_site,
            "init" => init,
            "term" => term,
            "elong" => elong,
            "km" => km,
            "kp" => kp,
            "kp_list" => kp_list,
            "Conc" => conc, 
            "deg_t" => deg_t,
            "deg_g" => deg_g,
            "run_time" => run_time,
            "starting_t" => starting_t,
            "delta_t" => delta_t
        ),
        "Results" => Dict(item => Dict() for item in outputs)
        )
end
dataL850_tau120_paper_log = save_data_to_dict(
        outputs, 
        saving_dict, 
        results850_tau120_paper_log,
        kp_list)

@save("dataL850_tau120_paper_log.jld", dataL850_tau120_paper_log)




@save("dataL350_tau120_paper.jld", dataL350_tau120_paper)
@save("dataL100_tau120_paper.jld", dataL100_tau120_paper)
@save("dataL850_tau120_paper.jld", dataL850_tau120_paper)
@save("dataL1350_tau120_paper.jld", dataL1350_tau120_paper)



L, l_ribosome, track_site = 100, 10, 4
ell = l_ribosome
ϵ = 20.0
init, term, elong = 0.2, ϵ, ϵ
deg_t = deg_g = 1/120
conc = range(0.1, stop = 12, length = 24)
km = 0.084 / 60
kp =  0.034 / 60 
kp_list = collect(kp .* conc)


kp_list = 10 .^ collect(-5:0.1:-1)

results_paper_time_pos_entry_L100_log = pmap(
    kp -> gather_first_pause_distribution(100, l_ribosome, track_site, init, term, elong, 1/120,km, kp, 100_000),
    kp_list)
@save("results_paper_time_pos_entry_L100_log.jld", results_paper_time_pos_entry_L100_log)

pos_test = [v[1] for v in results_paper_time_pos_entry_L100_log]

Plots.scatter(log10.(kp_list), pos_test ./100, ylims= (0,1), xlabel="log10(kp)", ylabel="Average position of first pause (s)", title="Average position of first pause vs kp", legend=false)

@time results_paper_time_pos_entry_L350_log = pmap(
    kp -> gather_first_pause_distribution(350, l_ribosome, track_site, init, term, elong, 1/120,km, kp, 20_000 * 3),
    kp_list)
@save("results_paper_time_pos_entry_L350_log.jld", results_paper_time_pos_entry_L350_log)

results_paper_time_pos_entry_L850_log = pmap(
    kp -> gather_first_pause_distribution(850, l_ribosome, track_site, init, term, elong, 1/120,km, kp, 20_000 * 4),
    kp_list)
@save("results_paper_time_pos_entry_L850_log.jld", results_paper_time_pos_entry_L850_log)

results_paper_time_pos_entry_L1350_log = pmap(
    kp -> gather_first_pause_distribution(1350, l_ribosome, track_site, init, term, elong, 1/120,km, kp, 20_000 * 5),
    kp_list)
@save("results_paper_time_pos_entry_L1350_log.jld", results_paper_time_pos_entry_L1350_log)
