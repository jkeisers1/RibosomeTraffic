using Random
using Distributions
using Distributed
using JLD
using ProfileView
using Interact
using Plots
using Blink

num_cores = Sys.CPU_THREADS

# Add workers (e.g., num_cores - 1)
addprocs(max(1, num_cores - 1))

@everywhere include("Gillespie_obc.jl")
@everywhere include("Basic_functions_obc.jl")
@everywhere include("analyitcal_functions.jl")
include("density_profile.jl")


surivial(tp, R, kp) = exp(-0.5*R*kp*tp^2)
surivial_single(kp,elong,init, L) = exp(-kp*first_passage_time(L, init, elong))
first_passage_time(L, init, elong) = L/(elong - init)

@everywhere begin
    L, l_ribosome, track_site = 350, 10, 1

    ϵ = 20.0
    init, term, elong = 0.2, 10*ϵ, ϵ
    deg_t = deg_g = 0

    k₋ = 0.084 / 60    
    k₊ =  0.034 / 60 

    kp = k₊
    km = k₋
end


km
######################################### different Lenght set kp and initiation #################################################################
L_list = collect(50:100:1450)

kp = 10^-3

@time P_L_list5 = pmap(L->gather_hit_probability(
    L, l_ribosome, track_site, init, term, elong, 
    km, kp, 10_000), L_list)

s_single_L3 = surivial_single.(kp,elong,init, L_list)
s_mulitple_L3 = surivial.(first_passage_time.(L_list, init, elong), J_LD_ext(init, elong, l_ribosome), kp)

plot(
    ylims = (0,1.05), 
    yticks = 0:0.2:1,
    guidefont = 22,
    xlims = (-75,1500),
    xticks = (0:300:1500),          # Font size for x and y labels (guides)
    tickfont = 18,           # Font size for tick labels
    legendfont = 18,
    ylabel = "S(tfp)",
    xlabel = "mRNA length",
    xrotation = 45,
    legend = :bottomleft
    )

scatter!(L_list,P_L_list5, markersize = 8, markercolor = colors[1], markershape = markers[1], label = "Numerics")
plot!(L_list, s_single_L3, linewidth = 6, linecolor = colors[end-2], label = "Literature")
plot!(L_list, s_mulitple_L3, linewidth = 6, linecolor = colors[3], label= "Model")
ylims!(0,1.05)

savefig("surivial_fraction_vs_gene_length_kp0001_init025.pdf")
    