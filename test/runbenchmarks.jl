# Pkg deps
using BenchmarkTools, SimpleDifferentialOperators, Test, Statistics

# Benchmark groups
benchmarks = BenchmarkGroup()

# Put in specific benchmarks
benchmarks["uniformReflecting"] = @benchmarkable diffusionoperators($(range(0.0, 1.0, length = 500)), $(Reflecting()), $(Reflecting()))
benchmarks["uniformMixed"] = @benchmarkable diffusionoperators($(range(0.0, 1.0, length = 500)), $(Mixed(ξ = 0.)), $(Mixed(ξ = 1.)))
benchmarks["irregularReflecting"] = @benchmarkable diffusionoperators($(collect(range(0.0, 1.0, length = 500))), $(Reflecting()), $(Reflecting()))
benchmarks["uniformMixed"] = @benchmarkable diffusionoperators($(collect(range(0.0, 1.0, length = 500))), $(Mixed(ξ = 0.)), $(Mixed(ξ = 1.)))

# Run and condense benchmarks
results = run(benchmarks)
results = median(results)

# To save results, manually call in the REPL: BenchmarkTools.save("benchmarks.json", results)

#Compare to old results
try
  oldresults= BenchmarkTools.load("benchmarks.json")[1]
  judge(oldresults, results)
catch err
  error("Couldn't load file- make sure that you've previously saved results.", err.prefix)
end
