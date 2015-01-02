using Compat, BucklingOfPipes

forcesAndMoments =
 @compat Dict(
    :Nx => 300, 
    :Ny => 100,
    :Nxy => 200,
    :Mx => 5, 
    :My => 3,
    :Mxy => 1)

materialProperties =@compat Dict(
  :E1 => 20.0*10^6, 
  :E2 => 1.5*10^6,
  :G12 => 1.0*10^6,
  :ν12 => .34,
  :alpha1 => 0.2e-6,
  :alpha2 => 20e-6,
  :beta1 => 0.2e-4,
  :beta2 => 20e-4)

laminateProperties =@compat Dict(
  :tk => 0.0075,
  :nplies => 5,
  :repeats => 1,
  :symmetric => false,
  :orientation => [0, 45, 90, 45, -45],
  :deltaTemp => -270.0,
  :deltaMoisture => 101
)

println()
createLaminate!(laminateProperties, materialProperties, forcesAndMoments)
show(laminateProperties)

println("ABD matrix:")
round(laminateProperties[:ABD], 3) |> display

println()
println("Total forces and moments:")
round(laminateProperties[:fTotal], 7) |> display

