using Random
using Distributions
using Distributed
using StatsBase
using JLD
using Pkg
using Plots
using ForwardDiff
using Pkg
using CairoMakie
using Colors
@everywhere include("Gillespie_obc.jl")
@everywhere include("Basic_functions_obc.jl")
@everywhere include("analyitcal_functions.jl")

L, l_ribosome, track_site = 300, 10, 4

ϵ = 16
init, term, elong = 0.02 , ϵ, ϵ
deg_t = deg_g = 0


k₋ = 0.084 / 60
k₊ = 0.034 / 60  * 1
#kp = k₊
#KD = 0.5
#k₋ = KD * k₊
#RNAP pausing
#k₋ = 1#/ ϵ
#k₊ = 0.1 #/ ϵ
#=
#tet
KD = 5 #μM the range given in papers is 0.5 - 20 μM
k₊ = 0.285 / ϵ # μM per second measured by tritton et.al 
k₋ = KD * k₊ / ϵ
=#
lowest_rate = min(k₋, k₊, round(mean(init),digits=3), elong, term)
    
run_time, starting_t, delta_t = 1*10^6, 10^5, 500

time, kymo, intern = Gillespie_obc_kymo_time(L, l_ribosome, init, term, elong, k₋, k₊, run_time, starting_t, delta_t)

fig = plot_kymo2(kymo', intern',k₊, 1/k₋, time, l_ribosome, track_site)


ylims!(0,50000)

fig = plot_kymo2_heatmap(kymo', intern', k₊, 1/k₋, time, l_ribosome, track_site;
                         fig_size=(2000,1000), yflip=true, rasterize_body=true)

p = plot_kymo2(kymo', intern', k₊, 1/k₋, time, l_ribosome, track_site;
                fig_size=(2000,1200), guidefontsize=20, tickfontsize=12)

xticks!(0:100:300)
xlims!(0,300)


savefig("kymo_translation_L$(L)_kp$(k₊)_km$(k₋)_init$(round(init, digits=4))_kymo.png")

kymo

sum(kymo'[101,1:end])

time

intern
test = (L * (L-1) / 2) * (kp/(elong+kp))^2 *(elong/(elong + kp))^(L-2)

savefig("kymo_translation_L$(L)_kp$(k₊)_km$(k₋)_init$(round(init, digits=4)).png")
