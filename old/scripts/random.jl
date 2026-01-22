
@save "data_L1350_tau300.jld" data_L1350_tau300 
@save "data_100_all.jld" data_100_all
@save "data_1350_all.jld" data_1350_all

@load "data_100_all.jld" data_100_all
@load "data_1350_all.jld" data_1350_all


################ T0 ##############################
Pos_100 = [v[1] for v in results_time_pos_100]
time_100 = [v[2] for v in results_time_pos_100]

T0_1350 = calculate_T0.(init, 20, 1350, kp_list, 120)
T0_100 = calculate_T0.(init, 20, 100, kp_list, 120)

Plots.plot(log10.(kp_list), T0_1350, label = "L: 1350", lw = 6)
Plots.plot!(log10.(kp_list), T0_100, label = "L: 100", lw = 6)
Plots.scatter!(log10.(kp_list), time_1350, label = false, markersize = 6, xlabel = "log(kp)", ylabel = "Time until first pause T0")
Plots.scatter!(log10.(kp_list), time_100, label = false, markersize = 6)
Plots.ylims!(0, 150)
savefig("T0_vs_kp.svg")


############# Tp ##############################
Tp_1350 = calculate_Tp_total.(1350, kp_list, km, l_ribosome, init, elong, 120)

Tp_100 = calculate_Tp_total.(100, kp_list, km, l_ribosome, init, elong, 120)

Tp_1350_test = 120 .- T0_1350


Plots.scatter!(log10.(kp_list), Tp_1350)
Plots.scatter!(log10.(kp_list), Tp_1350_test)

Plots.plot(log10.(kp_list), Tp_1350, label = "L: 1350", lw = 6)
Plots.plot!(log10.(kp_list), Tp_100, label = "L: 100", lw = 6)
Plots.scatter!(log10.(kp_list), time_1350, label = "Time until first pause", lw = 6, xlabel = "log(kp)", ylabel = "Time (s)")
Plots.scatter!(log10.(kp_list), time_100, label = false, lw = 6)
Plots.ylims!(0, 150)
savefig("Tp_vs_kp.png")


sort(data_1350_all["Results"]["density"])
paused_frac = [v[2]/v[4] for (k,v) in sort(data_1350_all["Results"]["density"])]
jammed_frac = []

data_100_all["Input"]

P0_100_num = [v[1]/v[3] for (k,v) in sort(data_100_all["Results"]["total_time"])]
P0_1350_num = [v[1]/v[3] for (k,v) in sort(data_1350_all["Results"]["total_time"])]
Pp_100_num = [v[2]/v[3] for (k,v) in sort(data_100_all["Results"]["total_time"])]
Pp_1350_num = [v[2]/v[3] for (k,v) in sort(data_1350_all["Results"]["total_time"])]

J100_num = [v[3] for (k,v) in sort(data_100_all["Results"]["current"])]
J1350_num = [v[3] for (k,v) in sort(data_1350_all["Results"]["current"])]


densities_1350 = calculate_state_densities_corrected_plateau4.(1350, kp_list, km, l_ribosome, init, elong, 120)
densities_100 = calculate_state_densities_corrected_plateau4.(100, kp_list, km, l_ribosome, init, elong, 120)


densities_1350 = calculate_state_densities_corrected_plateau.(1350, kp_list, km, l_ribosome, init, elong, 120)
densities_100 = calculate_state_densities_corrected_plateau.(100, kp_list, km, l_ribosome, init, elong, 120)

du100 = calc_rho_unpaused_contribution.(init, l_ribosome, elong, 100, kp_list, 120)
dp100 = calc_density_paused_total.(km, kp_list, l_ribosome, 100, init, elong, 120)
du1350 = calc_rho_unpaused_contribution.(init, l_ribosome, elong, 1350, kp_list, 120)
dp1350 = calc_density_paused_total.(km, kp_list, l_ribosome, 1350, init, elong, 120)

P0_100 = [v[1] for v in densities_100]
d0_100 = [v[2] for v in densities_100]
Pp_100 = [v[3] for v in densities_100]
dp_100 = [v[4] for v in densities_100]
dtot_100 = P0_100 .* d0_100 .+ Pp_100 .* dp_100

P0_1350 = [v[1] for v in densities_1350]
d0_1350 = [v[2] for v in densities_1350]
Pp_1350 = [v[3] for v in densities_1350]
dp_1350 = [v[4] for v in densities_1350]
dtot_1350 = P0_1350 .* d0_1350 .+ Pp_1350 .* dp_1350

#dp100_ana = calc_rho_unpaused_contribution.(init, l_ribosome, elong, 100, kp_list, 120)
#du100_ana = calc_density_paused_total.(km, kp_list, l_ribosome, 100, init, elong, 120)
#d100_ana = calculate_density_final_strict.(init, l_ribosome, elong, 100, kp_list, km, 120)

#du1350_ana = calc_rho_unpaused_contribution.(init, l_ribosome, elong, 1350, kp_list, 120)
#dp1350_ana = calc_density_paused_total.(km, kp_list, l_ribosome, 1350, init, elong, 120)
#d1350_ana = calculate_density_final_strict.(init, l_ribosome, elong, 1350, kp_list, km, 120)

Plots.scatter(ylabel = "P0", xlabel = "log(kp)", legend = :topright)
Plots.scatter!(log10.(kp_list), P0_100_num, xticks= -5:1: -1, label = "L:100, Paused Numerical", markersize = 6)
Plots.scatter!(log10.(kp_list), Pp_100_num, xticks= -5:1: -1, label = "L:100, Unpaused Numerical", markersize = 6)
Plots.plot!(log10.(kp_list), Pp_100, xticks= -5:1: -1, label = "L:100, Paused Analytical", lw = 6)
Plots.plot!(log10.(kp_list), P0_100, xticks= -5:1: -1, label = "L:100, Unpaused Analytical", lw = 6)

Plots.plot(log10.(kp_list), Pp_1350, xticks= -5:1: -1, label = "Paused Ana", lw = 6)
Plots.plot!(log10.(kp_list), P0_1350, xticks= -5:1: -1, label = "Unpaused Ana", lw = 6)
Plots.scatter!(log10.(kp_list), P0_1350_num, xticks= -5:1: -1, label = "unpaused num", markersize = 6)
Plots.scatter!(ylims = (0,1.05),yticks=0:0.5:1)
Plots.scatter!(log10.(kp_list), Pp_1350_num, xticks= -5:1: -1, label = "paused num", lw = 6)
savefig("P0_vs_kp_total.svg")


Plots.plot(log10.(kp_list), P0_100 .* d0_100  , xticks= -5:1: -1, label = "Paused Ana", lw = 6)
Plots.plot!(log10.(kp_list), Pp_100 .* dp_100 , xticks= -5:1: -1, label = "Unpaused Ana", lw = 6)
Plots.scatter!(log10.(kp_list), du100_num, xticks= -5:1: -1, label = "unpaused num", lw = 6)
Plots.scatter!(log10.(kp_list), dp100_num, xticks= -5:1: -1, label = "paused num", lw = 6)

Plots.scatter(ylabel = "Density", xlabel = "log(kp)", ylims = (0,0.042), yticks=0:0.02:0.04)
Plots.scatter!(log10.(kp_list), d1350_num, xticks= -5:1: -1, label = "L:1350, Paused Numerical", markersize = 6)
Plots.plot!(log10.(kp_list), dtot_1350, xticks= -5:1: -1, label = "L:1350, total Analytical", lw = 6)
Plots.scatter!(log10.(kp_list), d100_num, xticks= -5:1: -1, label = "L:100, Paused Numerical", markersize = 6)
Plots.plot!(log10.(kp_list), dtot_100, xticks= -5:1: -1, label = "L:100, total Analytical", lw = 6)
savefig("density_vs_kp_total.svg")



d100_num = [v[4] for (k,v) in sort(data_100_all["Results"]["density"])]
dp100_num = [v[6] for (k,v) in sort(data_100_all["Results"]["density"])]
du100_num = [v[5] for (k,v) in sort(data_100_all["Results"]["density"])]

d1350_num = [v[4] for (k,v) in sort(data_1350_all["Results"]["density"])]
dp1350_num = [v[6] for (k,v) in sort(data_1350_all["Results"]["density"])]
du1350_num = [v[5] for (k,v) in sort(data_1350_all["Results"]["density"])]

first_approx100 = weight_tau_cutoff.(100, kp_list, k₋, l_ribosome, init, elong, 120)
first_approx1350 = weight_tau_cutoff.(1350, kp_list, k₋, l_ribosome, init, elong, 120)

Plots.plot(conc_ana, P0_100, lw = 6, label= "L100, new analytical")
Plots.plot!(conc_ana, T0_ana1350 ./ (T0_ana1350 .+ Tp_ana1350), lw = 6, label = false)
Plots.scatter!(ylabel = "P0", xlabel = "Conc", legend = :topright)
Plots.scatter!(log10.(kp_list), P_100, ylims = (0,1.05), yticks = 0:0.5:1, xticks= -5:1: -1, label = "L = 100")
#Plots.plot!(conc_ana, first_approx100, lw = 4, ls = :dash, label = "L100, First Approximation")

Plots.scatter!(log10.(kp_list), P_1350, ylims = (0,1), yticks = 0:0.5:1, xticks= -5:1: -1, label = "L = 1350")
Plots.plot!(conc_ana, T0_ana1350 ./ (T0_ana1350 .+ Tp_ana1350), lw = 6, label = "L1350, new analytical")
#Plots.plot!(conc_ana, first_approx1350, lw = 4, ls = :dash, label = "L1350, First Approximation")
savefig("P0_vs_conc_zoomed.png")

Tp_ana100 = calculate_Tp_total.(100, kp_list, k₋, l_ribosome, init, elong, 120)
T0_ana100 = calculate_T0.(init, elong, 100, kp_list, 120)

Tp_ana1350 = calculate_Tp_total.(1350, kp_list, k₋, l_ribosome, init, elong, 120)
T0_ana1350 = calculate_T0.(init, elong, 1350, kp_list, 120)

Tp_ana100 .+ T0_ana100

Pos_100_ana = calculate_avg_pause_position_robust.(init, elong, 100, kp_list, 120)

Pos_1350_ana = calculate_avg_pause_position_robust.(init, elong, 1350, kp_list, 120)

kp_list

Plots.plot(ylabel = "Average pausing position", xlabel = "log(kp)", legend = :topright, ylims = (0,1), yticks=0:0.5:1)
Plots.plot!(log10.(kp_list), Pos_100_ana ./ 100, xticks= -5:1: -1, label = "L = 100", lw = 6, ylims = (0,1))
Plots.scatter!(log10.(kp_list), Pos_100 ./ 100, xticks= -5:1: -1, label = "L = 100",  markersize = 6)

    
Plots.plot!(log10.(kp_list), Pos_1350_ana ./ 1350, xticks= -5:1: -1, label = "L = 1350", lw = 6)
Plots.scatter!(log10.(kp_list), Pos_1350 ./ 1350, xticks= -5:1: -1, label = "L = 1350", markersize = 6)
savefig("avg_pausing_position_vs_kp.svg")


states_100 = calculate_state_currents.(100, kp_list, k₋, l_ribosome, init, elong, 120)
states_1350 = calculate_state_currents.(1350, kp_list, k₋, l_ribosome, init, elong, 120)

P0_100 = [v[1] for v in states_100]
J0_100 = [v[2] for v in states_100]
Pp_100 = [v[3] for v in states_100]
Jp_100 = [v[4] for v in states_100]
Jtot_100 = P0_100 .* J0_100 .+ Pp_100 .* Jp_100

P0_1350 = [v[1] for v in states_1350]
J0_1350 = [v[2] for v in states_1350]
Pp_1350 = [v[3] for v in states_1350]
Jp_1350 = [v[4] for v in states_1350]
Jtot_1350 = P0_1350 .* J0_1350 .+ Pp_1350 .* Jp_1350

Plots.scatter(ylabel = "P0", xlabel = "log(kp)", legend = :topright, ylims = (0,1))
Plots.plot!(log10.(kp_list), T0_100 ./ (T0_100 .+ Tp_100), xticks= -5:1: -1, label = "L:100, Analytical", lw = 6)
Plots.plot!(log10.(kp_list), P0_1350, xticks= -5:1: -1, label = "L:1350, Analytical", lw = 6)
Plots.scatter!(log10.(kp_list), P0_100_num, xticks= -5:1: -1, label = "L:100, Numerical", lw = 6)
Plots.scatter!(log10.(kp_list), P0_1350_num, xticks= -5:1: -1, label = "L:1350, Numerical", lw = 6)
savefig("P0_vs_kp.png")

init
Plots.scatter(ylabel = "J", xlabel = "log(kp)", legend = :topright, ylims = (0,0.21), yticks=0:0.05:0.2)
Plots.plot!(log10.(kp_list), Jtot_100, label = "L:100, Analytical", linewidth = 8)
Plots.plot!(log10.(kp_list), Jtot_1350, label = "L:1350, Analytical", linewidth = 8)
#Plots.plot!(log10.(kp_list), P0_1350 .* J0_1350 .+ Jp1350_2, xticks= -5:1: -1, label = "L:1350, Analytical (Instant Fill)", lw = 6)
Plots.scatter!(log10.(kp_list), J100_num, label = "L:100, Numerical", markersize = 6)
Plots.scatter!(log10.(kp_list), J1350_num, label = "L:1350, Numerical", markersize = 6)
Plots.ylims!(-0.01,0.21)
savefig("J_vs_kp.svg")

states_100 = calculate_state_densities.(100, kp_list, k₋, l_ribosome, init, elong, 120)
states_350 = calculate_state_densities.(350, kp_list, k₋, l_ribosome, init, elong, 120)
states_1350 = calculate_state_densities.(1350, kp_list, k₋, l_ribosome, init, elong, 120)

P0_100 = [v[1] for v in states_100]
d0_100 = [v[2] for v in states_100]
Pp_100 = [v[3] for v in states_100]
dp_100 = [v[4] for v in states_100]
dtot_100 = P0_100 .* d0_100 .+ dp100


P0_1350 = [v[1] for v in states_1350]
d0_1350 = [v[2] for v in states_1350]
Pp_1350 = [v[3] for v in states_1350]
dp_1350 = [v[4] for v in states_1350]
dtot_1350 = P0_1350 .* d0_1350 .+ dp1350_2
dtot_1350_2 = P0_1350 .* d0_1350 .+ dp1350_ins

P0_350 = [v[1] for v in states_350]
d0_350 = [v[2] for v in states_350]
Pp_350 = [v[3] for v in states_350]
dp_350 = [v[4] for v in states_350]
dtot_350 = P0_350 .* d0_350 .+ dp350


#dp100 = calc_density_jam_only_exact.(km, kp_list, l_ribosome, 100, init, elong, 120)
dp100 = calc_rho_paused_final.(km, kp_list, l_ribosome, 100, init, elong, 120)
dp1350_ins = calc_density_jam_only_instant_fill.(km, kp_list, l_ribosome, 1350, init, elong, 120)

dp2 = calculate_density_jam_only_exact2.(km, kp_list, l_ribosome, 100, init, elong, 120)

dp100 .- dp2

dp350 = calc_rho_paused_final.(km, kp_list, l_ribosome, 350, init, elong, 120)
dp1350_2 = calc_rho_paused_final.(km, kp_list, l_ribosome, 1350, init, elong, 120)
dp1350 = calc_density_jam_only_exact.(km, kp_list, l_ribosome, 1350, init, elong, 120)


Plots.scatter(ylabel = "Density", xlabel = "log(kp)")
Plots.plot!(conc_test[18:29], dtot_100[18:29], label = "L:100, paused Analytical", lw = 6)
Plots.plot!(conc_test[18:29], dtot_1350[18:29], label = "L:1350, paused Analytical", lw = 6)
Plots.plot!(conc_test[18:29], dtot_1350_2[18:29], label = "L:1350, paused Analytical (Instant Fill)", lw = 6)
Plots.plot!(conc_test[18:29], dtot_350[18:29], label = "L:350, paused Analytical", lw = 6)
Plots.scatter!(conc_test[18:29], d100_num[18:29], label = "L:100, Numerical", markersize = 6)
Plots.scatter!(conc_test[18:29], d1350_num[18:29], label = "L:1350, Numerical", markersize= 6)
Plots.ylims!(0,0.025)
savefig("density_vs_kp.png")

dens_ana = calculate_density_final.(init, l_ribosome, elong, 100, kp_list, k₋, 120)
dens_ana1350 = calculate_density_final.(init, l_ribosome, elong, 1350, kp_list, k₋, 120)


dens_p_100 = calc_rho_paused_corrected.(km, kp_list, l_ribosome, 100, init, elong, 120)
dens_p_1350 = calc_rho_paused_corrected.(km, kp_list, l_ribosome, 1350, init, elong, 120)

dens_u_100 = calc_rho_unpaused_part_exact.(init, l_ribosome, elong, 100, kp_list, 120)
dens_u_1350 = calc_rho_unpaused_part_exact.(init, l_ribosome, elong, 1350, kp_list, 120)


dens_u_100_num = [v[5] for (k,v) in sort(data_100_all["Results"]["density"])]
dens_p_100_num = [v[6] for (k,v) in sort(data_100_all["Results"]["density"])]

dens1350_num = [v[4] for (k,v) in sort(data_1350_all["Results"]["density"])]    

Plots.plot(log10.(kp_list), dens_p_100, xticks= -5:1: -1, label = "Analytical", lw = 6)
Plots.scatter!(log10.(kp_list), dens_p_100_num, xticks= -5:1: -1, label = "Simulation", lw = 6)

Plots.scatter(log10.(kp_list), dens_u_100, xticks= -5:1: -1, label = "Numerical", lw = 6)
Plots.plot!(log10.(kp_list), dens_ana, xticks= -5:1: -1, label = "Analytical L100", lw = 6)

Plots.plot(log10.(kp_list), dens_ana1350, xticks= -5:1: -1, label = "Analytical L1350", lw = 6)
Plots.scatter!(log10.(kp_list), dens1350_num, xticks= -5:1: -1, label = "Simulation L1350", lw = 6)

savefig("density_vs_kp.png")



Ju1350 = calc_J_unpaused_contribution.(init, elong, L, kp_list, l_ribosome, 120)
Jp1350 = calc_J_paused_part_exact.(km, kp_list, l_ribosome, 1350, init, elong, 120)
Ju1350_num = [v[1] for (k,v) in sort(data_1350_all["Results"]["current"])]
Jp1350_num = [v[2] for (k,v) in sort(data_1350_all["Results"]["current"])]
J1350_tot = [v[3] for (k,v) in sort(data_1350_all["Results"]["current"])]


Ju100 = calc_J_unpaused_part.(init, l_ribosome, elong, 100, kp_list, 120)
Jp100 = calc_J_paused_part_exact.(km, kp_list, l_ribosome, 100, init, elong, 120)
Ju100_num = [v[1] for (k,v) in sort(data_100_all["Results"]["current"])]
Jp100_num = [v[2] for (k,v) in sort(data_100_all["Results"]["current"])]
J100_tot = [v[3] for (k,v) in sort(data_100_all["Results"]["current"])]


Plots.scatter(ylabel = "J unpaused state", xlabel = "log(kp)", legend = :topright)
Plots.scatter(log10.(kp_list), Ju1350_num, xticks= -5:1: -1, label = "L:1350, Simulation", lw = 6, markercolor = :blue)
Plots.plot!(log10.(kp_list), Ju1350, xticks= -5:1: -1, label = "L=1350, Analytical", lw = 6, linecolor = :blue)

Plots.scatter(ylabel = "J paused state", xlabel = "log(kp)", legend = :topright)
Plots.scatter(log10.(kp_list), Jp1350_num, xticks= -5:1: -1, label = "L:1350, Simulation", lw = 6, markercolor = :blue)
Plots.plot!(log10.(kp_list), Jp1350, xticks= -5:1: -1, label = "L=1350, Analytical", lw = 6, linecolor = :blue)

Plots.scatter(ylabel = "J total", xlabel = "log(kp)", legend = :topright)
Plots.scatter!(log10.(kp_list), J1350_tot, xticks= -5:1: -1, label = "L:1350, Simulation", lw = 6, markercolor = :blue)
Plots.plot!(log10.(kp_list), Jp1350 .+ Ju1350, xticks= -5:1: -1, label = "L=1350, Analytical", lw = 6, linecolor = :blue)
Plots.scatter!(log10.(kp_list), J100_tot, xticks= -5:1: -1, label = "L:100, Simulation", lw = 6, markercolor = :blue)
Plots.plot!(log10.(kp_list), Jp100 .+ Ju100, xticks= -5:1: -1, label = "L=100, Analytical", lw = 6, linecolor = :blue)
Plots.xlims!(-4,-2)

savefig("current_vs_kp_zoom.png")

log10.(collect(k₊ .* conc))


Plots.scatter!(log10.(kp_list), J1350_num, xticks= -5:1: -1, label = "L:1350, Simulation", lw = 6, markercolor = :black)
Plots.plot!(log10.(kp_list), J1350, xticks= -5:1: -1, label = "L=1350,Analytical", lw = 6, linecolor = :black)


J100_num = [v[3] for (k,v) in sort(data_100_all["Results"]["current"])]

J100 = calculate_current_final.(init, l_ribosome, elong, 100, kp_list, k₋, 120)

J1350 = calculate_current_final.(init, l_ribosome, elong, 1350, kp_list, k₋, 120)



savefig("current_adjusted_vs_kp.png")

res = Gillespie_obc(L, l_ribosome, track_site, deg_t, deg_g, init, term, elong, k₋, kp_list[20], run_time, starting_t, delta_t)

tmp1 = res[1]
tmp2 = res[1]

deg_list = collect(range(0, stop = 1/60, length= 20))


run_Gillespie_for_concentration(outputs, x, Conc, k₋, k₊, L, l_ribosome, track_site, init, elong, term,deg_t, deg_g)

conc
results0 = run_Gillespie_for_concentration(outputs, 1000, conc, k₋, k₊, 400, l_ribosome, track_site, init, elong, term,0, 0)
results60 = run_Gillespie_for_concentration(outputs, 1000, conc, k₋, k₊, 400, l_ribosome, track_site, init, elong, term,0, 1/60)
results120 = run_Gillespie_for_concentration(outputs, 1000, conc, k₋, k₊, 400, l_ribosome, track_site, init, elong, term,0, 1/120)
results300 = run_Gillespie_for_concentration(outputs, 1000, conc, k₋, k₊, 400, l_ribosome, track_site, init, elong, term,0, 1/300)


P_0_0 = [v[2]/v[3] for (k,v) in sort(results0["Results"]["total_time"])]
P_0_60 = [v[2]/v[3] for (k,v) in sort(results60["Results"]["total_time"])]
P_0_120 = [v[2]/v[3] for (k,v) in sort(results120["Results"]["total_time"])]
P_0_300 = [v[2]/v[3] for (k,v) in sort(results300["Results"]["total_time"])]

Plots.scatter(conc, P_0_0, ylims = (0,1), yticks = 0:0.5:1, xticks= (0:2:12))
Plots.scatter!(conc, P_0_60, ylims = (0,1), yticks = 0:0.5:1, xticks= (0:2:12))
Plots.scatter!(conc, P_0_120, ylims = (0,1), yticks = 0:0.5:1, xticks= (0:2:12))

Plots.scatter!(conc, P_0_300, ylims = (0,1), yticks = 0:0.5:1, xticks= (0:2:12))


prot_syn0 = [v[3] for (k,v) in sort(results0["Results"]["current"])]
prot_syn60 = [v[3] for (k,v) in sort(results60["Results"]["current"])]
prot_syn120 = [v[3] for (k,v) in sort(results120["Results"]["current"])]
prot_syn300 = [v[3] for (k,v) in sort(results300["Results"]["current"])]


plot(conc, prot_syn0 ./ prot_syn0[1], lw = 6)
plot!(conc, prot_syn60 ./ prot_syn60[1], lw = 6)
plot!(conc, prot_syn120 ./ prot_syn120[1], lw =6)
plot!(conc, prot_syn300 ./ prot_syn300[1], lw = 6)
xticks!(0:2:12)
yticks!(0:0.5:1)
savefig("protein_syn_vs_pausing.svg")

plot!(conc[2:end], expected_proteins_test.(init, elong, 120, kp_list[2:end], L, l_ribosome))

results2 = run_Gillespie_for_deg_list(outputs, 5000, 2, k₋, k₊, L, l_ribosome, track_site, init, elong, term, deg_list)
results6 = run_Gillespie_for_deg_list(outputs, 500, 6, k₋, k₊, L, l_ribosome, track_site, init, elong, term, deg_list)
results10 = run_Gillespie_for_deg_list(outputs, 500, 10, k₋, k₊, L, l_ribosome, track_site, init, elong, term, deg_list)


xaxis = collect(keys(sort(results["Results"]["current"])))
prot_syn0 = [v[3] for (k,v) in sort(results0["Results"]["current"])]
prot_syn2 = [v[3] for (k,v) in sort(results2["Results"]["current"])]
prot_syn6 = [v[3] for (k,v) in sort(results6["Results"]["current"])]
prot_syn10 = [v[3] for (k,v) in sort(results10["Results"]["current"])]


plot(deg_list, prot_syn2)
plot!(deg_list, prot_syn6)
plot!(deg_list, prot_syn10)

ylims!(0,0.01)
scatter!(conc, prot_syn120)
scatter!(conc, prot_syn300)

plot(conc, savitzky_golay(prot_syn60 ./ prot_syn0, 5, 2).y ,label = "Fold change 60", lw = 6)
plot!(conc, savitzky_golay(prot_syn120 ./ prot_syn0, 5, 2).y ,label = "Fold change 120", lw = 6)
plot!(conc, savitzky_golay(prot_syn300 ./ prot_syn0, 5, 2).y,label = "Fold change 300", lw =6)

ylims!(0,30)
xticks!(0:2:12)

savefig("required_fold_change.svg")


function movmean_simple(data, window)
    n = length(data)
    half = floor(Int, window ÷ 2)
    result = zeros(Float64, n)
    for i in 1:n
        left = max(1, i - half)
        right = min(n, i + half)
        result[i] = mean(data[left:right])
    end
    return result
end

function backward_movmean(data, w)
    n = length(data)
    result = zeros(Float64, n)
    for i in 1:n
        left = max(1, i-w+1)
        result[i] = mean(data[left:i])
    end
    return result
end



function run_Gillespie_for_init_list(outputs, x, α_list, k₋, k₊, L, l_ribosome, track_site, elong, term)
    
    lowest_rate = min(k₋, k₊, round(mean(α_list),digits=3), elong, term)
    
    tss = round((1/lowest_rate)) * 100
    delta_t = round((1/lowest_rate))*x

    #tss = round((1/10)) *x
    #delta_t = round((1/10))*x*10

    starting_t, delta_t = tss, delta_t
    run_time = starting_t + delta_t
    @time results = pmap(
    α -> Gillespie_obc(L, l_ribosome, track_site, deg_t, deg_g, α, term, elong, k₋, k₊, run_time, starting_t, delta_t; kymo=false),α_list)

    saving_dict = Dict(
    "Input" => Dict(
        "L" => L,
        "l_ribosome" => l_ribosome,
        "track_site" => track_site,
        "α_list" => α_list,
        "term" => term,
        "elong" => elong,
        "k₋" => k₋,
        "k₊" => k₊, 
        "deg_t" => deg_t,
        "deg_g" => deg_g,
        "run_time" => run_time,
        "starting_t" => starting_t,
        "delta_t" => delta_t
    ),
    "Results" => Dict(item => Dict() for item in outputs)
    )

    data = save_data_to_dict(
        outputs, 
        saving_dict, 
        results, 
        α_list
    )

    return data
end

function run_Gillespie_for_concentration(outputs, x, Conc, k₋, k₊, L, l_ribosome, track_site, init, elong, term,deg_t, deg_g)
    
    k₊_list = Conc .* k₊
    tss = round((1/k₊)) 
    delta_t = round((1/k₊))* x

    #tss = round((1/10)) *x
    #delta_t = round((1/10))*x*10

    starting_t, delta_t = tss, delta_t
    run_time = starting_t + delta_t

    @time results = pmap(
    kp -> Gillespie_obc(L, l_ribosome, track_site, deg_t, deg_g, init, term, elong, k₋, kp, run_time, starting_t, delta_t; kymo=false),
    k₊_list)

    saving_dict = Dict(
    "Input" => Dict(
        "L" => L,
        "l_ribosome" => l_ribosome,
        "track_site" => track_site,
        "init" => init,
        "term" => term,
        "elong" => elong,
        "k₋" => k₋,
        "k₊" => k₊,
        "Conc" => Conc, 
        "deg_t" => deg_t,
        "deg_g" => deg_g,
        "run_time" => run_time,
        "starting_t" => starting_t,
        "delta_t" => delta_t
    ),
    "Results" => Dict(item => Dict() for item in outputs)
    )

    data = save_data_to_dict(
        outputs, 
        saving_dict, 
        results, 
        k₊_list
    )

    return data
end



function run_Gillespie_for_deg_list(outputs, x, Conc, k₋, k₊, L, l_ribosome, track_site, init, elong, term, deg_list)
    
    kp = Conc * k₊
    tss = round((1/k₊)) 
    delta_t = round((1/k₊))* x

    #tss = round((1/10)) *x
    #delta_t = round((1/10))*x*10

    starting_t, delta_t = tss, delta_t
    run_time = starting_t + delta_t

    @time results = pmap(
    deg -> Gillespie_obc(L, l_ribosome, track_site, 0, deg, init, term, elong, k₋, kp, run_time, starting_t, delta_t; kymo=false),
    deg_list)

    saving_dict = Dict(
    "Input" => Dict(
        "L" => L,
        "l_ribosome" => l_ribosome,
        "track_site" => track_site,
        "init" => init,
        "term" => term,
        "elong" => elong,
        "k₋" => k₋,
        "k₊" => k₊,
        "Conc" => Conc,
        "deg_list" => deg_list, 
        "deg_t" => deg_t,
        "deg_g" => deg_g,
        "run_time" => run_time,
        "starting_t" => starting_t,
        "delta_t" => delta_t
    ),
    "Results" => Dict(item => Dict() for item in outputs)
    )

    data = save_data_to_dict(
        outputs, 
        saving_dict, 
        results, 
        deg_list
    )

    return data
end


α_list = 10 .^ collect(-4:0.25:1)

α_list[4]

transcription_data0 = []
L_list0 = collect(60:60:600)


for L in L_list0
    t_data = run_Gillespie_for_init_list(outputs, 10000, α_list, k₋, k₊, L, l_ribosome, track_site, elong, term)
    push!(transcription_data0, t_data)
end

transcription_data1 = []
L_list1 = collect(700:100:1500)


for L in L_list1
    t_data = run_Gillespie_for_init_list(outputs, 8000, α_list, k₋, k₊, L, l_ribosome, track_site, elong, term)
    push!(transcription_data1, t_data)
end


transcription_data2 = []
L_list2 = collect(1800:200:3200)

for L in L_list2
    t_data = run_Gillespie_for_init_list(outputs, 6000, α_list, k₋, k₊, L, l_ribosome, track_site, elong, term)
    push!(transcription_data2, t_data)
end

transcription_data3 = []
L_list3 = collect(3400:200:5000)

for L in L_list3
    t_data = run_Gillespie_for_init_list(outputs, 6000, α_list, k₋, k₊, L, l_ribosome, track_site, elong, term)
    push!(transcription_data3, t_data)
end

p = scatter()
for item in transcription_data2
    tmp = extract_current_efficiency(item)
    p = scatter(log10.(α_list),tmp, ylims = (0,1))
    display(p)
end
p

tmp_tot = vcat(transcription_data0, transcription_data1, transcription_data2, transcription_data3)

transcription_data0

init070,length070 = obtain_current_eff(tmp_tot[2:end], 0.90)
init085,length085 = obtain_current_eff(tmp_tot[2:end], 0.95)
init0975,length0975 = obtain_current_eff(tmp_tot[2:end], 0.99)


scatter!(length070, init070, label = "0.7", legend=:topleft)
scatter!(length085, init085, label = "0.85")
ß


tmp = extract_current_efficiency(transcription_data1[end-2])
scatter(log10.(α_list),tmp)
ylabel!("Current efficiency")
xlabel!("log10(α)")
hline!([0.96])
str = replace("transcription_data_km$(k₋)_kp$(round(k₊,digits=5))", "." => "") * ".jld"
transcription_data = vcat(transcription_data0, transcription_data1, transcription_data2)
@save str transcription_data



#=
translation_data = Dict()

Conc_large = collect(range(0.5,stop=12, length = 11))
Conc_small = collect(range(0,stop=0.48, length = 9))
Conc = collect(0:0.5:12)


init = 0.002
pre_run = run_Gillespie_for_concentration(outputs, 0.001, Conc, k₋, k₊, L, l_ribosome, track_site, init, elong, term)
translation_data_kp1 = run_Gillespie_for_concentration(outputs, 50, Conc, k₋, k₊, L, l_ribosome, track_site, init, elong, term)

translation_data[init] = translation_data_kp1


init = 0.02
translation_data_kp2 = run_Gillespie_for_concentration(outputs, 25, Conc, k₋, k₊, L, l_ribosome, track_site, init, elong, term)

translation_data[init] = translation_data_kp2


init = 0.2
translation_data_kp3 = run_Gillespie_for_concentration(outputs, 10, Conc, k₋, k₊, L, l_ribosome, track_site, init, elong, term)
translation_data[init] = translation_data_kp3



init = 2
translation_data_kp4 = run_Gillespie_for_concentration(outputs, 5, Conc, k₋, k₊, L, l_ribosome, track_site, init, elong, term)
translation_data[init] = translation_data_kp4


str = replace("translation_data_L$(L)_km$(k₋)_kp$(round(k₊,digits=5))", "." => "") * ".jld"

@save str translation_data
=#