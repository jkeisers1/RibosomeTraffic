
Tfinal = 120.0  # or some max time
Δt = 1
bin_edges = collect(0:Δt:Tfinal)
n_bins = length(bin_edges) - 1
# Prepare an accumulator
accumulator = zeros(n_bins, 1001)  # suppose L=1000
n_trials = 1000
number_of_steps = 30_000
for trial_idx in 1:n_trials
    # call the combined function
    binned_occupancy = gather_density_profile_and_bin(
        L, l_ribosome, track_site, init, term, elong,
        km, kp,
        number_of_steps,
        bin_edges
    )
    # sum into accumulator
    accumulator .+= binned_occupancy
end


# Then get average across trials
accumulator ./= n_trials


density_profileL1000_init2_elong20 = accumulator    

bar(density_profileL1000_init2_elong20)
density_profileL1000_init2_elong20

@save "density_profileL1000_init02_elong20.jld" density_profileL1000_init02_elong20

@load "density_profileL1000_init02_elong20.jld"

w = Window()
ui = @manipulate for t in bin_edges[2:end]
    # find index
    row_index = findfirst(x-> x == t ,bin_edges[2:end])
    # Select the current row
    current_row = density_profileL1000_init2_elong20[row_index, :]
    analyitcal = instantaneousDensityProfile_ext(L,l_ribosome, init, elong, t) 
    # Plot the current row
    plot(
        current_row,
        xlabel = "lattice site",
        ylabel = "density",
        label = "Row $row_index",
        legend = :topright
    )

    plot!(
        analyitcal,
        label = false,
        linewidth = 5
    )
end


body!(w, ui)

@everywhere begin
    function bin_time_data(all_data::Vector{Matrix{Float64}}, Δt::Float64)
    # 1) Find the maximum time across all trials
    T_max = maximum([maximum(d[1:end, 2]) for d in all_data])  # d[:,2] is the time column

    # 2) Define time bins from 0 to T_max

    edges = 0:Δt:T_max
    n_bins = length(edges) - 1

    # 3) Figure out how many sites (L) from the first dataset
    #    (assuming all trials have the same shape)
    L = size(all_data[1], 2) - 2  # subtract 2 for (dt, t) columns
    n_trials = length(all_data)

    # We'll store:
    #   binned_occupancy[trial, bin, site]
    # and a bin_counts[trial, bin] for how many snapshots fell in that bin
    binned_occupancy = zeros(n_trials, n_bins, L)
    bin_counts       = zeros(n_trials, n_bins)

    # 4) Loop over trials
    for i in 1:n_trials
        data = all_data[i]
        # Each row in `data` is: [dt, t, occ1, occ2, ..., occL]
        times     = data[:, 2]
        occupancy = data[:, 3:end]  # shape: (N_steps x L)

        # For each time value, figure out which bin index it belongs to
        # We'll use `searchsortedlast` to find bin edges
        bin_indices = searchsortedlast.(Ref(edges), times)
        # This can return values up to n_bins+1, clamp them to [1, n_bins]
        bin_indices = clamp.(bin_indices, 1, n_bins)

        # Accumulate occupancies in that bin
        for (row_idx, b) in enumerate(bin_indices)
            binned_occupancy[i, b, :] .+= occupancy[row_idx, :]
            bin_counts[i, b] += 1
        end
    end

    bin_counts
    print("hello")
    # 5) Convert accumulated sums to means (snapshot-average) within each bin
    for i in 1:n_trials
        for b in 1:n_bins
            if bin_counts[i, b] > 0
                binned_occupancy[i, b, :] ./= bin_counts[i, b]
            end
        end
    end

    binned_occupancy

    # 6) Now average over the trials
    # binned_occupancy has shape: (n_trials x n_bins x L)
    final_occupancy = mean(binned_occupancy, dims=1)  # result: 1 x n_bins x L
    final_occupancy
    final_occupancy = dropdims(final_occupancy, dims=1)  # shape: (n_bins x L)

    # We often want the bin centers (rather than edges) for plotting
    bin_centers = (edges[1:end-1] .+ edges[2:end]) ./ 2

    return bin_centers, final_occupancy
    end

    function bin_single_run(
    profile::Matrix{Float64}, 
    bin_edges::Vector{Float64})
    # profile columns: [dt, t, occ1, occ2, ..., occL]
    times     = profile[:, 2]          # shape: (Nsteps,)
    occupancy = profile[:, 3:end]      # shape: (Nsteps, L)

    n_bins = length(bin_edges) - 1
    L      = size(occupancy, 2)
    binned_matrix = zeros(n_bins, L)
    bin_counts    = zeros(n_bins)

    # For each snapshot time, figure out which bin it belongs to
    # We'll use `searchsortedlast` to find the appropriate bin index
    bin_indices = searchsortedlast.(Ref(bin_edges), times)
    # clamp them to [1, n_bins] in case time == last edge, etc.
    bin_indices = clamp.(bin_indices, 1, n_bins)

    # Accumulate occupancy
    for (i, b) in enumerate(bin_indices)
        @inbounds binned_matrix[b, :] .+= occupancy[i, :]
        bin_counts[b] += 1
    end

    # Convert sums to average occupancy per bin
    for b in 1:n_bins
        if bin_counts[b] > 0
            binned_matrix[b, :] ./= bin_counts[b]
        end
    end

    return binned_matrix
    end

    function gather_density_profile_and_bin(
    L::Int, l_ribosome::Int, track_site, init::Float64, term::Float64, elong::Float64,
    k_minus::Float64, k_plus::Float64,
    number_of_steps::Int64,
    bin_edges::Vector{Float64})
    # 1) Get the single-run profile (large matrix) from Gillespie
    profile = gather_density_profile(
        L, l_ribosome, track_site, init, term, elong,
        km, kp, number_of_steps
    )

    # 2) Immediately bin it
    binned_result = bin_single_run(profile, bin_edges)

    # 3) Return ONLY the binned result (not the huge profile)
    return binned_result
    end

    function run_Gillespie_for_init_list(outputs, x, α_list, k₋, k₊, L, l_ribosome, track_site, elong, term)
    
    lowest_rate = min(k₋, k₊, round(mean(α_list),digits=3), elong, term)
    
    tss = round((1/lowest_rate)) *x
    delta_t = round((1/lowest_rate))*x*10

    #tss = round((1/10)) *x
    #delta_t = round((1/10))*x*10

    starting_t, delta_t = tss, delta_t
    run_time = starting_t + delta_t

    @time results = pmap(
    α -> Gillespie_obc(L, l_ribosome, track_site, deg_t, deg_g, α, term, elong, k₋, k₊, run_time, starting_t, delta_t; kymo=false),
    α_list)

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

    function run_Gillespie_for_concentration(outputs, x, Conc, k₋, k₊, L, l_ribosome, track_site, init, elong, term)
    
    k₊_list = Conc .* k₊
    tss = round((1/k₊)) *100
    delta_t = round((1/k₊))*1000 * x

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
end