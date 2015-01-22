using Compat, ClassicalLaminateTheory

forcesAndMoments =
 @compat Dict(
    :Nx => 1000, 
    :Ny => 300,
    :Nxy => 150,
    :Mx => 2.5, 
    :My => 5,
    :Mxy => 10)

materialProperties =@compat Dict(
  :E1 => [20.0e6, 13.0*10^6],
  :E2 => [1.5*10^6, 1.0*10^6],
  :G12 => [1.0*10^6, 0.75*10^6],
  :ν12 => [0.34, 0.34],
  :alpha1 => [0.2e-6, 0.2e-5],
  :alpha2 => [0.2e-4, 0.2e-3],
  :beta1 => [0.2e-6, 0.2e-5],
  :beta2 => [0.2e-4, 0.2e-3],
  :thickness => [0.0075, 0.0055])

laminateProperties =@compat Dict(
  :nplies => 4,
  :materials => [2, 2, 1, 1],
  :repeats => 2,
  :symmetric => true,
  :symmetricrepeats => 1,
  :orientation => [0, 90, 0, 90],
  :deltaTemp => -280.0,
  :deltaMoisture => 102
)

qm = qmat(materialProperties)
materialProperties |> display

println("\nQ:")
qm |> display

qb0 = qbarmat(qm, 0.0)
qb90 = qbarmat(qm, 90.0)
println("\nQbar 0 Degrees:")
qb0 |> display
println("\nQbar 90 Degrees:")
qb90 |> display

println()
createLaminate!(laminateProperties, materialProperties, forcesAndMoments)
show(laminateProperties)

println("ABD matrix:")
round(laminateProperties[:ABD], 3) |> display

println()
println("Total forces and moments:")
round(laminateProperties[:fTotal], 7) |> display

