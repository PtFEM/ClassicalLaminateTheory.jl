# Top level test script for Stan.jl
using ClassicalLaminateTheory
using Base.Test

println("Running tests for ClassicalLaminateTheory-j0.3/4-v0.0.1:")

tests = [
  "test_clt1.jl",
  "test_clt2.jl",
  "test_clt3.jl",
  "test_clt4.jl"
]

println("Running tests:")

for my_test in tests
    println("\n  * $(my_test) *")
    include(my_test)
end

println()