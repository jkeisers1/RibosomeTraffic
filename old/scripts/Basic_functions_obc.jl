#using Random
#using Distributions
#using Pkg
#using Distributed
#using BenchmarkTools

function create_lattice_l(L::Int64, l_ribosome::Int64) 

    #init_region = zeros(Int64, l_ribosome+1)
    # lattice = zeros(Int64, L)
    #lattice = zeros(Int64, L-2*l_ribosome)

    #end_region = zeros(Int64, l_ribosome)
    
    #lattice = vcat(init_region, lattice, end_region)

    return zeros(Int64, L+1)
end


# function create_lattice(L, ρ)
#     lattice = zeros(L+1)
#     while mean(lattice) < ρ
#         index = findfirst(x -> x == 0, lattice)
#         lattice[index] += 1
#     end
    
#     shuffle!(lattice)
#     lattice[1] = 0
#     return lattice
# end
function create_rates_l(init, term, elong, L)  #strength of bottleneck is set to uniform rate by default
       
    #rates = zeros(1+l_ribosome+L) .+ elong # for checking
    rates = zeros(1+L) .+ elong
    rates[1]= init   #setting initation rate to α at site 1
    
    rates[end] = term #setting termination rate to β at site L+1
    
    return rates
end

# function create_rates(L, ϵ; α = 1.0, β = 1.0)
#     rates = ones(L+1) .* ϵ
#     rates[1] = α
#     rates[end] = β
#     return rates
# end

# L, l_ribosome= 1, 10
# init, term, elong = 1, 1, 1
# bn = Dict()
# lattice = create_lattice_l(L, l_ribosome)
# rates = create_rates_l(init, term, elong, bn, L, l_ribosome)
# posR = positionRibosomes(lattice)::Vector{Int64}

function positionRibosomes(lattice)::Vector{Int64}
    position = []
    for j in 1:length(lattice)
        if lattice[j] == 1 
            push!(position, j)
        else
            push!(position, 0)
        end
    end
    position[1] = 1
    return position
end
# L, ρ = 10, 0.5
# init, elong, term = 100, 1, 0.2
# lattice =  create_lattice(L, ρ)
# rates = create_rates(L, elong, α=init, β=term)
# posR = positionRibosomes(lattice)
# elongation_vector = get_elongation_vector_obc(lattice, rates)
# hcat(lattice, posR, elongation_vector)

# moving_particle = sample(posR, Weights(elongation_vector))

function create_log10_range_list(starting_point, stopping_point, step_expo)
    list = []
    tmp = collect(range(starting_point,stop=stopping_point, step=step_expo))
    for item in tmp
        push!(list,10^(item))
    end
    return list
end

function create_log2_range_list(starting_point, stopping_point, step_expo)
    list = []
    tmp = collect(range(starting_point,stop=stopping_point, step=step_expo))
    for item in tmp
        push!(list,2^(item))
    end
    return list
end
