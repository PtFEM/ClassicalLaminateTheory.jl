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
  :E2 => 1.5*10^6,
  :G12 => 1.0*10^6,
  :ν12 => .34,
  :alpha1 => 0.2e-6,
  :alpha2 => 0.2e-4,
  :beta1 => 0.2e-6,
  :beta2 => 0.2e-4,
  :thickness => 0.0075)

laminateProperties =@compat Dict(
  :nplies => 10,
  :repeats => 1,
  :symmetric => true,
  :orientation => [90, 45, -45, 90, 45, 45, 90, -45, 45, 90],
  :deltaTemp => 0.0,
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

