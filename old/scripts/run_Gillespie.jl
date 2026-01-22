using Random
using Distributions
using Distributed
using Plots
using JLD
using Colors
number_of_workers = 4
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


function times_from_dict(dict)
    tw = collect(values(sort(Dict( key => value[1] for (key,value) in sort(dict["Results"]["total_time"])))))
    tc = collect(values(sort(Dict( key => value[2] for (key,value) in sort(dict["Results"]["total_time"])))))
    t = collect(values(sort(Dict( key => value[3] for (key,value) in sort(dict["Results"]["total_time"])))))
    return tw, tc, t
end


function densities_from_dict(dict)
    m = collect(values(sort(Dict( key => value[1] for (key,value) in sort(dict["Results"]["density"])))))
    p = collect(values(sort(Dict( key => value[2] for (key,value) in sort(dict["Results"]["density"])))))
    j = collect(values(sort(Dict( key => value[3] for (key,value) in sort(dict["Results"]["density"])))))
    tot = collect(values(sort(Dict( key => value[4] for (key,value) in sort(dict["Results"]["density"])))))
    tot_unp = collect(values(sort(Dict( key => value[5] for (key,value) in sort(dict["Results"]["density"])))))
    tot_p = collect(values(sort(Dict( key => value[6] for (key,value) in sort(dict["Results"]["density"])))))
    return m, p, j, tot, tot_unp, tot_p
end


function current_from_dict(dict)
    Jw = collect(values(sort(Dict( key => value[1] for (key,value) in sort(dict["Results"]["current"])))))
    Jc = collect(values(sort(Dict( key => value[2] for (key,value) in sort(dict["Results"]["current"])))))
    J = collect(values(sort(Dict( key => value[3] for (key,value) in sort(dict["Results"]["current"])))))
    return Jw, Jc, J
end


function cluster_length_from_dict(dict)
    return sort(dict["Results"]["cluster_length_distribution"])
end

function cluster_number_unp_from_dict(dict)
    return sort(dict["Results"]["cluster_number_distribution_unp"])
end

function cluster_number_p_from_dict(dict)
    return sort(dict["Results"]["cluster_number_distribution_p"])
end

function paused_part_distr_from_dict(dict)
    return sort(dict["Results"]["paused_part_distribution"])
end

L, l_ribosome, track_site = 1000, 10, 4
init, term, elong = 0.2, 20, 1

k₊_list = 10.0 .^ collect(range(-4, stop = -1, length = 6))
deg_t = deg_g = 0
k₋ = 1.4 * 10^-3
k₊ = 1

lowest_rate = 1/0.6
starting_t, delta_t = lowest_rate *100, lowest_rate * 1000
run_time = starting_t + delta_t

x = 4
run_time, starting_t, delta_t = 12 *10^x, 4*10^x, 8*10^x 
@time c1 = Gillespie_obc(L, l_ribosome, track_site, deg_t, deg_g, init, term, elong, k₋, k₊, run_time, starting_t, delta_t)


#define outputs
outputs = [
    "current", 
    "density", 
    "paused_part_distribution", 
    "cluster_length_distribution", 
    "cluster_number_distribution_unp",
    "cluster_number_distribution_p",
    "total_time",
    ]


# set simulation data for saving
β_crit

x = 5
run_time, starting_t, delta_t = 20 * 10^x, 5*10^x, 15*10^x


@time results_LD500_k₋001_k₊0001 = pmap(
    α -> Gillespie_obc(L, l_ribosome, track_site, deg_t, deg_g, α, term, elong, k₋, k₊, run_time, starting_t, delta_t; kymo=false),
    α_list)

simulation_data_LD500_k₋001_k₊0001 = Dict(
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

data_LD500_k₋001_k₊0001 = save_data_to_dict(outputs, simulation_data_LD500_k₋001_k₊0001, results_LD500_k₋001_k₊0001,α_list)

@save "data_LD500_k₋001_k₊0001.jld" data_LD500_k₋001_k₊0001

dist_num = cluster_number_from_dict(data_LD500_k₋001_k₊0001)


bar(dist_num[α_list[40]], xlims = (-1,100))

@save "data_LD250_k₋01_k₊001.jld" data_LD250_k₋01_k₊001

@load "data_LD250_k₋01_k₊001.jld"
@load "data_HD_250_k₋1_k₊1.jld"
@load "data_LD_250_k₋1_k₊1.jld"

dens_LD2 = density_from_dict(data_LD)
current_LD2 = current_from_dict(data_LD)
dens_LD = density_from_dict( data_LD500_k₋001_k₊0001)
current_LD = current_from_dict(data_LD500_k₋001_k₊0001)

scatter(dens_LD[1:30], current_LD[1:30])
vline!([ρ_max_fct(k₋,k₊,1)])

@time results_HD250_k₋01_k₊001 = pmap(
    β -> Gillespie_obc(L, l_ribosome, track_site, deg_t, deg_g, init, β, elong, k₋, k₊, run_time, starting_t, delta_t; kymo=false),
    β_list)

simulation_data_HD250_k₋01_k₊001 = Dict(
    "Input" => Dict(
        "L" => L,
        "l_ribosome" => l_ribosome,
        "track_site" => track_site,
        "init" => init,
        "β_list" => β_list,
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

data_HD250_k₋01_k₊001 = save_data_to_dict(outputs, simulation_data_HD250_k₋01_k₊001, results_HD250_k₋01_k₊001, β_list)



dens_HD = density_from_dict( data_HD250_k₋01_k₊001)
current_HD = current_from_dict(data_HD250_k₋01_k₊001)
dens_HD2 = density_from_dict(data_HD)
current_HD2 = current_from_dict(data_HD)

scatter(dens_HD, current_HD)
scatter!(dens_LD, current_LD)
scatter!(dens_HD2, current_HD2)
scatter!(dens_LD2, current_LD2)
@save "data_HD250_k₋01_k₊001.jld" data_HD250_k₋01_k₊001

ρL_1000_HD = density_from_dict(data_HD_L1000_k₋0001_k₊001)
JL_1000_HD = current_from_dict(data_HD_L1000_k₋0001_k₊001)

scatter(ρL_1000_HD,JL_1000_HD )

ρL_500_HD = density_from_dict(data_HD_L500_k₋0001_k₊001)
JL_500_HD = current_from_dict(data_HD_L500_k₋0001_k₊001)

# set simulation data for saving
simulation_data_LD = Dict(
    "Input" => Dict(
        "L" => L,
        "l_ribosome" => l_ribosome,
        "track_site" => track_site,
        "init" => init,
        "α_list" => α_list,
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


@time results_LD_L1000_k₋0001_k₊001 = pmap(
    α -> Gillespie_obc(L, l_ribosome, track_site, deg_t, deg_g, α, term, elong, k₋, k₊, run_time, starting_t, delta_t; kymo=false),
    α_list)

data_LD_L1000_k₋0001_k₊001 = save_data_to_dict(outputs, simulation_data_LD, results_LD_L1000_k₋0001_k₊001, α_list)

@save "data_LD_L1000_k₋0001_k₊001.jld" data_LD_L1000_k₋0001_k₊001
@save "data_HD_L1000_k₋0001_k₊001.jld" data_HD_L1000_k₋0001_k₊001

@load "data_LD_L1000_k₋0001_k₊001.jld"
@load "data_HD_L1000_k₋0001_k₊001.jld"

@load "data_βcrit_Llist.jld"

@load "data_L150.jld"
@load "data_L100.jld"
@load "data_L50.jld"

data_βcrit_Llist["Results"]

L_list = data_βcrit_Llist["Input"]["L_list"]
dens_β = density_from_dict(data_βcrit_Llist)

scatter(
    L_list, 
    dens_β,
    ylabel = "density",
    xlabel = "System size",
    label = "β = β_crit",
    markersize = 5, 
    ylims = (0.2,0.32)
    )
hline!([ρ_max_fct(0.01,0.001,1)],linewidth = 4, label = "ρ max analytical")


β_num_dist= cluster_number_from_dict(data_βcrit_Llist)


savefig("density_vs_system_size_β_crit.pdf")

bar(
    number_distribution[1], 
    xlabel = "Number of clusters", 
    ylabel = "Distribution",
    label = "$(k₊_list[1])",
    xlims = (0,40)
    )

savefig("bar_k+-6.pdf")

bar(
    number_distribution[15], 
    xlabel = "Number of clusters", 
    ylabel = "Distribution",
    label = "k₊= $(round(k₊_list[15],digits=5))",
    xlims = (0,40)
    )

savefig("bar_k+-4.pdf")

bar(
    number_distribution[21], 
    xlabel = "Number of clusters", 
    ylabel = "Distribution",
    label = "k₊= $(round(k₊_list[21],digits=5))",
    xlims = (0,30)
    )

    savefig("bar_k+-3.pdf")    
bar(
    number_distribution[end], 
    xlabel = "Number of clusters", 
    ylabel = "Distribution",
    label = "k₊= $(round(k₊_list[end],digits=5))",
    xlims = (0,30)
    )
    savefig("bar_k+-1.pdf") 
number_distribution = data_L150["cluster_number_distribution"]
k₊_list
scatter(k₊_list, yscale =:log10)
k₊_list = data_L150["k₊_list"]
k₊_list_axis = log10.(k₊_list)
data_L150["k₊_list"][1:16]


ρ150 = [item[4] for item in data_L150["density"]]
ρ100 = [item[4] for item in data_L100["density"]]
ρ50 =  [item[4] for item in data_L50["density"]]

J150 = [item[1] for item in data_L150["current"]]
J100 = [item[1] for item in data_L100["current"]]
J50 =  [item[1] for item in data_L50["current"]]

colors = [
    RGB(249/255, 65/255, 68/255),
    RGB(243/255, 114/255, 44/255),
    RGB(248/255, 150/255, 30/255),
    RGB(249/255, 132/255, 74/255),
    RGB(249/255, 199/255, 79/255),
    RGB(144/255, 190/255, 109/255),
    RGB(67/255, 170/255, 139/255),
    RGB(77/255, 144/255, 142/255),
    RGB(87/255, 117/255, 144/255),
    RGB(39/255, 125/255, 161/255)
]

marker_shapes = [:circle, :star5, :diamond, :square, :utriangle]
line_styles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
label_font_size = 20
legend_font_size = 12
axis_font_size = 18
tick_label_fontsize = 18
markersize = 6
linewidth = 5
Plots.plot(
    #title = "Sample Plot",
    xlabel = "k₊ (log)\n",
    ylabel = "density",
    ylims = (0, 0.6),
    yticks = (0:0.2:0.6),
    titlefontsize = label_font_size,
    legendfontsize = legend_font_size,
    xguidefontsize = axis_font_size,  # Font size for x-axis label
    yguidefontsize = axis_font_size,
    xtickfontsize = tick_label_fontsize,
    ytickfontsize = tick_label_fontsize,
    #xlims = (0,1),
    #xticks = (0:0.2:1)
    )

scatter!(k₊_list_axis, ρ_max_fct.(10^-3, k₊_list,1), label = "ρ_max pTASEP",markersize = markersize, markercolor = colors[end], markershape=marker_shapes[4])

scatter!(k₊_list_axis, ρ50, label = "L = 50", markersize = markersize, markercolor = colors[1], markershape=marker_shapes[1])
scatter!(k₊_list_axis, ρ100, label ="L = 100",markersize = markersize,markercolor = colors[4], markershape=marker_shapes[2])
scatter!(k₊_list_axis, ρ150, label ="L = 150",markersize = markersize,markercolor = colors[5], markershape=marker_shapes[3])

plot!(k₊_list_axis, open_bio_light.(0.1,k₊_list, 10^-3, 50,1), label = false ,linewidth= linewidth,linecolor = colors[1])
plot!(k₊_list_axis, open_bio_light.(0.1,k₊_list, 10^-3, 100,1),label = false ,linewidth= linewidth,linecolor = colors[2])
plot!(k₊_list_axis, open_bio_light.(0.1,k₊_list, 10^-3, 150,1),label = false ,linewidth= linewidth, linecolor = colors[5])


savefig("ρ_vs_k₊_OBC.pdf")

scatter

scatter

dens = density_from_dict(data_βcrit_Llist)
L_list = data_βcrit_Llist["Input"]["L_list"]
scatter(L_list,dens)
hline!([0.25])
data_LD_L1000_k₋0001_k₊001["Results"]["cluster_number_distribution"]
ρL_1000_LD = density_from_dict(data_LD_L1000_k₋0001_k₊001)
JL_1000_LD = current_from_dict(data_LD_L1000_k₋0001_k₊001)

β_list = data_HD_L1000_k₋0001_k₊001["Input"]["β_list"]
length_distr = cluster_length_from_dict(data_HD_L1000_k₋0001_k₊001)
number_distr = cluster_number_from_dict(data_HD_L1000_k₋0001_k₊001)
ρL_1000_HD = density_from_dict(data_HD_L1000_k₋0001_k₊001)
JL_1000_HD = current_from_dict(data_HD_L1000_k₋0001_k₊001)
scatter(ρL_1000_HD,JL_1000_HD)
scatter(length_distr[β_list[1]],xlims =(0,500))
findmin(ρL_1000_HD)

scatter(ρL_1000_LD,JL_1000_LD)
scatter(ρL_1000_HD,JL_1000_HD)
scatter!(ρL_500_LD,JL_1000_LD)
scatter!(ρL_500_HD,JL_1000_HD)
vline!([ρ_max_fct(k₋, k₊,1)])
ρL_500_LD = density_from_dict(data_LD_L500_k₋0001_k₊001)
JL_500_LD = current_from_dict(data_LD_L500_k₋0001_k₊001)

ρL_250_LD = density_from_dict(data_LD_L250)
JL_250_LD = current_from_dict(data_LD_L250)

ρL_250_HD = density_from_dict(data_HD_L250)
JL_250_HD = current_from_dict(data_HD_L250)

sort(data_LD_L500_k₋0001_k₊001["Results"]["density"])

scatter(ρL_500_LD,JL_500_LD)

function density_from_dict(dict)
    return collect(values(sort(Dict( key => value[4] for (key,value) in sort(dict["Results"]["density"])))))
end
function current_from_dict(dict)
    return collect(values(sort(dict["Results"]["current"])))
end
function cluster_length_from_dict(dict)
    return sort(dict["Results"]["cluster_length_distribution"])
end
function cluster_number_from_dict(dict)
    return sort(dict["Results"]["cluster_number_distribution"])
end


#@save "data_HD_250_k₋001_k₊0001.jld" data_HD_L250
#@save "data_LD_250_k₋001_k₊0001.jld" data_LD_L250

@load "data_HD_250_k₋1_k₊1.jld"
@load "data_LD_250_k₋1_k₊1.jld"

@load "data_HD_250_k₋001_k₊0001.jld"
@load "data_LD_250_k₋001_k₊0001.jld"

colors = palette(:roma25)

dens_HD_250 = collect(values(sort(Dict( key => value[4] for (key,value) in sort(data_HD_L250["Results"]["density"])))))
dens_LD_250 = collect(values(sort(Dict( key => value[4] for (key,value) in sort(data_LD_L250["Results"]["density"])))))


current_HD_250 = collect(values(sort(data_HD_L250["Results"]["current"])))
current_LD_250 = collect(values(sort(data_LD_L250["Results"]["current"])))

α_crit_fct(k₋, k₊,1)

α_crit = ρ_max_fct(k₋, k₊,1) * fa_hat(k₋, k₊,ρ_max_fct(k₋, k₊,1) , 1)

α_list_LD = data_LD_L250["Input"]["α_list"]
β_list_LD = data_HD_L250["Input"]["β_list"]
scalefontsizes(1.5)
scatter(
    β_list,
    JL_500_HD,
    label = " OBC numerics", 
    xlabel = "β", 
    ylabel = "density ",
    title = "k₋ = $(round(k₋,digits =2)),k₊ = $(round(k₊ ,digits =3))"
    )

scatter(
    α_list,
    ρL_500_LD,
    label = " OBC numerics", 
    xlabel = "α", 
    ylabel = "density",
    title = "k₋ = $(round(k₋,digits =2)),k₊ = $(round(k₊ ,digits =3))"
)

hline!([ρ_max_fct(k₋, k₊, 1)], label ="ρ_max")

vline!([α_crit], label = "α crit", linewidth = 5)

vline!([1 - β_crit], label = "β crit from MC to HD", linewidth = 5)
hline!([Wang(ρ_max_fct(k₋, k₊, 1), k₋, k₊,1)], label = "ρ_max", linewidth = 5)
savefig("critical_values_β.pdf")
ρ_max_fct(k₋, k₊, 1)
ρ_max_fct(k₋, k₊, 1)*(1-fp(k₋, k₊))

label_font_size = 20
legend_font_size = 12
axis_font_size = 25
tick_label_fontsize = 18
Plots.plot(
    #title = "Sample Plot",
    xlabel = "density",
    ylabel = "current",
    #label = "Data Series 1",
    titlefontsize = label_font_size,
    legendfontsize = legend_font_size,
    xguidefontsize = axis_font_size,  # Font size for x-axis label
    yguidefontsize = axis_font_size,
    xtickfontsize = tick_label_fontsize,
    ytickfontsize = tick_label_fontsize,
    xlims = (0,1.05),
    xticks = (0:0.2:1)
    )

k₋,k₊  = 10^-2, 10^-3 
α_list_ana = collect(range(0.001, stop = α_crit, length = 200))
β_list_ana = collect(range(0.001, stop = β_crit, length = 200))


plot!(
    ρ_LD.(α_list_ana, k₋, k₊, 1), 
    Wang_LD.(α_list_ana, k₋, k₊, 1), 
    linewidth = 4, 
    linecolor = colors[end],
    label = "analytical LD",
    )


plot!(
    ρ_HD.(β_list_ana, 1),
    Wang_HD.(β_list_ana, k₋, k₊, 1),
    linewidth = 4, 
    linecolor = colors[1],
    label = "analytical HD "
    )

scatter!(ρL_500_HD, JL_500_HD, label = "numerical HD", markersize = 5, markercolor = colors[1])
scatter!(ρL_500_LD[1:30], JL_500_LD[1:30], label = "numerical LD",markersize = 5, markercolor = colors[end])
scatter!(ρL_500_LD[1:30], JL_500_LD[1:30], label = "numerical LD",markersize = 5, markercolor = colors[end])
scatter!(ρL_1000_HD,JL_1000_HD )

vline!([ρ_max_fct(k₋, k₊, 1)], label = "ρ_max",linewidth =4)
savefig("obc_LD_HD.pdf")
ρL_500_LD
JL_500_LD
scatter(α_list[1:25], ρL_500_LD[1:25])

α_list_ana = collect(range(0.001, stop = α_crit_fct(1, 1, 1) , length = 200))
β_list_ana = collect(range(0.001, stop = β_crit_fct(1, 1, 1) , length = 200))
plot!(
    ρ_LD.(α_list_ana, 1, 1, 1), 
    Wang_LD.(α_list_ana, 1, 1, 1), 
    linewidth = 4, 
    linecolor = colors[5],
    linestyle = :dot,
    label = false,
    )

plot!(
    ρ_HD.(β_list_ana, 1),
    Wang_HD.(β_list_ana, 1, 1, 1),
    linewidth = 4, 
    linecolor = colors[7],
    linestyle = :dot,
    label = false
    )

scatter!(
    x_values_LD, 
    y_values_LD, 
    label = false,
    markersize = 5, 
    markercolor = colors[6],
    markershape = :star
    )
scatter!(
    x_values_HD, 
    y_values_HD, 
    label = false,
    markersize = 5, 
    markershape = :star,
    markercolor = colors[8])



savefig("obc_J_vs_dens_shifted_and_rescaled.pdf")

xlims!(0,1)
ylims!(0,0.006)
plot!(ρ_list_ana, Wang.(ρ_list_ana, k₋, k₊, 1), linewidth = 4)

ρ_max = ρ_max_fct(k₋,k₊, 1)

vline!([ρ_max])

@load "J_data_250.jld"
J250 = J 

@load "J_data_500.jld"
J500 = J

k₊_t, k₋_t = 0.01, 0.1

α_crit = α_crit_fct(k₋_t, k₊_t, 1)
β_crit = β_crit_fct(k₋_t, k₊_t, 1)

α_list_ana = collect(range(0.001, stop = α_crit , length = 200))
β_list_ana = collect(range(0.001, stop = β_crit, length = 200))


k₊_list = [10^-6, 10^-5, 10^-4, 10^-3, 0.01, 10^-1, 10^0, 10]
k₋_list = [10^-6, 10^-5, 10^-4, 10^-3, 0.01, 10^-1, 10^0, 10]

index2 = findfirst(x-> x==k₋_t, k₋_list)
index1 = findfirst(x-> x==k₊_t, k₊_list)

Plots.plot(
    #title = "Sample Plot",
    xlabel = "ρ",
    ylabel = "J",
    #label = "Data Series 1",
    titlefontsize = label_font_size,
    legendfontsize = legend_font_size,
    xguidefontsize = axis_font_size,  # Font size for x-axis label
    yguidefontsize = axis_font_size,
    xtickfontsize = tick_label_fontsize,
    ytickfontsize = tick_label_fontsize,
    xlims = (0,1.05),
    xticks = (0:0.2:1)
    )

scatter(
    ρ_list, J500[index2,index1][1:end,1], 
    label = "PBC Numerics", markersize = 5, markercolor = colors[15]
    )

vline!([ρ_max_fct(k₋_t, k₊_t,1) ], linewidth = 5)
scatter!(ρL_500_HD, JL_500_HD, label = "numerical HD", markersize = 5, markercolor = colors[1])
scatter!(ρL_500_LD[1:30], JL_500_LD[1:30], label = "numerical LD",markersize = 5, markercolor = colors[end])
scatter!(ρL_500_LD[1:30], JL_500_LD[1:30], label = "numerical LD",markersize = 5, markercolor = colors[end])
scatter!(ρL_250_HD, JL_250_HD)
scatter!(ρL_250_LD, JL_250_LD)
savefig("tmp.pdf")
    
scatter(β_list, ρL_500_HD, label = "L=250", ylims=(0,1))
vline!([β_crit*fa_hat(k₋, k₊, ρ_max_fct(k₋, k₊,1), 1) / fa(k₋, k₊)])
fa_hat(k₋, k₊, ρ_max_fct(k₋, k₊,1), 1) / fa(k₋, k₊)

β_crit*fa_hat(k₋, k₊, ρ_max_fct(k₋, k₊,1), 1) / fa(k₋, k₊)
scatter!(β_list, ρL_500_HD, label = "L=500")

index2 = findfirst(x-> x== k₋, k₋_list)
index1 = findfirst(x-> x== k₋, k₊_list)

scatter(
    ρ_list, J500[index2,index1][1:end,1], 
    label = false, markersize = 5, markercolor = colors[15], markershape = :star
    )

plot!(
    ρ_list, 
    Wang.(ρ_list, k₋_t, k₊_t, 1), 
    linewidth = 4, 
    linecolor = colors[3],
    label = false
    )
plot!(
    ρ_list, 
    Wang.(ρ_list, 1, 1, 1), 
    linewidth = 4, 
    linecolor = colors[3],
    linestyle =:dot,
    label = "Analytical"
    )



findmax(J250[index2,index1][1:end,1])
ρ_list[9]
ρ_max_fct(k₋_t, k₊_t,1)
savefig("pbc_J_vs_dens_shifted_and_rescaled.pdf")