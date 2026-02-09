using RibosomeTraffic
using BenchmarkTools
using JumpProcesses
using Printf

function run_benchmark()
    println("Starting Benchmark... (Please wait, silencing noisy outputs)")

    # 1. SETUP & WARMUP (The Noisy Part)
    # We create a "Black Hole" for text. Nothing can escape this block.
    # ---------------------------------------------------------------
    L = 200
    model = TranscriptModel(L, 10, 0.5, 1.0, ones(L), 0.01, 0.1, 0.0)
    t_max = 500.0
    t_sciml = 0.0
    t_custom = 0.0

    # Redirect stdout to null (the trash can)
    original_stdout = stdout
    redirect_stdout(devnull)

    try
        # Run warmup (this usually prints the massive text)
        run_sciml_simulation(model, 1.0)
        run_custom_simulation(model, 1.0)

        # Run benchmarks (BenchmarkTools prints progress, we hide that too)
        t_sciml = @belapsed run_sciml_simulation($model, $t_max)
        t_custom = @belapsed run_custom_simulation($model, $t_max)
    finally
        # ALWAYS restore standard output, even if there is an error
        redirect_stdout(original_stdout)
    end
    # ---------------------------------------------------------------

    # 2. PRINT RESULTS (Now that stdout is back)
    println("\n=========================================")
    println("   PERFORMANCE AUDIT: CUSTOM VS SCIML    ")
    println("=========================================")
    @printf "1. SciML Solver Time:  %.8f seconds\n" t_sciml
    @printf "2. Custom Solver Time: %.8f seconds\n" t_custom
    println("-----------------------------------------")

    speedup = t_sciml / t_custom
    @printf "   SPEEDUP FACTOR:    %.8fx faster\n" speedup
    println("=========================================")
end

run_benchmark()
