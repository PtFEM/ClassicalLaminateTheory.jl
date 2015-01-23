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
  :orientation => [0, 90, 0, 90],
  :materials => [2, 2, 1, 1],
  :repeats => 2,
  :symmetric => true,
  :symmetricrepeats => 1,
  :deltaTemp => -280.0,
  :deltaMoisture => 102
)

createLaminate!(laminateProperties, materialProperties, forcesAndMoments)

fm = [-402.962, -1102.96, 150.0, 2.5, 5.0, 10.0]

@assert(round(laminateProperties[:fTotal], 2) == round(fm, 2))
println("Tests for CLT3 passed!")
