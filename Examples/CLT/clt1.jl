using Compat, BucklingOfPipes

forcesAndMoments =
  @compat Dict(
    :Nx => 300, 
    :Ny => 100,
    :Nxy => 0,
    :Mx => 5, 
    :My => 0,
    :Mxy => 0) 

materialProperties = @compat Dict(
  :E1 => 20.0*10^6, 
  :E2 => 1.5*10^6,
  :G12 => 1.0*10^6,
  :ν12 => .34,
  :alpha1 => 0.2e-6,
  :alpha2 => 20e-6,
  :beta1 => 0.2e-4,
  :beta2 => 20e-4)

laminateProperties = @compat Dict(
  :tk => 0.0075,
  :nplies => 4,
  :repeats => 1,
  :symmetric => true,
  :orientation => [0, 45, 90, -45],
  :deltaTemp => 0.0,
  :deltaMoisture => 0.00
)

qm = qmat(materialProperties)

println("\nQ:")
qm |> display

qb0 = qbarmat(qm, 0.0)
qb90 = qbarmat(qm, 90.0)
qb45 = qbarmat(qm, 45.0)
qbm45 = qbarmat(qm, -45.0)
println("\nQbar 0 Degrees:")
qb0 |> display
println("\nQbar 90 Degrees:")
qb90 |> display
println("\nQbar 45 Degrees:")
qb45 |> display
println("\nQbar -45 Degrees:")
qbm45 |> display

println()
createLaminate!(laminateProperties, materialProperties, forcesAndMoments)
show(laminateProperties)

println("ABD matrix:")
round(laminateProperties[:ABD], 3) |> display

println()
println("Total forces and moments:")
round(laminateProperties[:fTotal], 7) |> display

