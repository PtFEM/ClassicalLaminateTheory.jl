using Compat, ClassicalLaminateTheory

forcesAndMoments =
 @compat Dict(
    :Nx => 0, 
    :Ny => 0,
    :Nxy => 0,
    :Mx => 0, 
    :My => 0,
    :Mxy => 0)

materialProperties =@compat Dict(
  :E1 => 20.0*10^6, 
  :E2 => 1.1*10^6,
  :G12 => 0.65*10^6,
  :ν12 => .21,
  :alpha1 => 0.2e-6,
  :alpha2 => 0.2e-4,
  :beta1 => 0.2e-6,
  :beta2 => 0.2e-4,
  :thickness => 0.0075)

laminateProperties =@compat Dict(
  :nplies => 4,
  :repeats => 2,
  :symmetric => true,
  :orientation => [0, 45, 90, -45],
  :deltaTemp => -280.0,
  :deltaMoisture => 0
)

println()
createLaminate!(laminateProperties, materialProperties, forcesAndMoments)
show(laminateProperties)

println("ABD matrix:")
round(laminateProperties[:ABD], 3) |> display

println()
println("Total forces and moments:")
round(laminateProperties[:fTotal], 7) |> display

