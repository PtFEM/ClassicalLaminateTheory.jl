# Top level test script for Stan.jl
using ClassicalLaminateTheory
using Base.Test

println("Running tests for ClassicalLaminateTheory-j0.3/4-v0.0.1:")

tests = [
  "test_clt.jl"
]

println("Running tests:")

for my_test in tests
    println("\n  * $(my_test) *")
    include(my_test)
end

println()