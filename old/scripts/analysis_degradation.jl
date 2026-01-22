using JLD
using Plots
using Colors


# Example usage:
num_colors = 10
colors = generate_distinguishable_colors(num_colors)
println(colors)

### loading data

@load "deg_J.jld" 
@load "deg_ρ.jld" 
@load "deg_P_res.jld"
@load "deg_paused_part.jld"
@load "deg_time.jld" 


# define functions 

ρ_profile(init, elong, pos, deg) = init / elong * exp(-deg * pos / velo(init, elong)) # density profile. The profile is assuming steady state density times an exponential decay function
avg_ρ(init, elong, v, δ, L) = init/elong * v / δ  * (1 - exp(-δ * L / v)) / L# integration of the density profile over the entire lattice
velo(init, elong) = elong * (1 - init/elong) # average speed of the particle on the lattice
J_deg(ρ, δ, v, β) = ρ * exp(-δ*L/v) * β # first part of the equation is the density times the degradation and exponential decay function but for i = L where i is the lattice site
Pw(δ, k₋, ρ, L, k₊) = (δ + k₋) / (δ + k₋ + ρ * L * k₊)

J_total(δ, k₋, k₊, β, v, ρ, L) = Pw(δ, k₋, ρ, L, k₊) * J_deg(ρ, δ, v, β)

# parameters used for simualtion --> add them to the data loading package

avg_ρ(init, elong, velo(init, elong), deg_list[2], L)


L, l_ribosome, track_site = 250, 1, 1
init, term, elong = 0.3, 10, 10
k₋, k₊ = 0.084 / 60 , 0.034 / 60

mRNA_life_time = reverse([1, 2, 3, 4, 5]) .* 60
deg_list = 1 ./ mRNA_life_time
pushfirst!(deg_list, 0 )

Conc = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12]
k₊_list = Conc .* k₊
deg_list
# plotting
#scalefontsizes(1.2)
colors = distinguishable_colors(6,[RGB(1,1,1), RGB(0,0,0)], dropseed = true)
scalefontsizes(1.2)
p = plot(ylabel = "J", xlabel = "Conc")
for i in 2:1:6
    δ = deg_list[i]
    tmp_ρ = sum(ρ_profile.(init, elong, 1:0.1:250, δ)) / length(1:0.1:250)
    analytical_current = J_total.(δ, k₋, k₊_list, term, velo(init, elong),init/elong, L)
    numerical_current = deg_J[i]
    conc_labels = string.(Conc)
    deg_labels = string.(round.(deg_list,digits = 4))
    label = deg_labels[i]
    scatter!(Conc, numerical_current, label = "δ: $label", markercolor = colors[i], markershape = :auto,markersize = 5)
    plot!(Conc, analytical_current, label = false, linecolor = colors[i], linewidth = 3)
end

p
ylims!(0,0.31)

savefig("J_vs_Conc.png")