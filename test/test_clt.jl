using Compat, ClassicalLaminateTheory

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
  :beta2 => 20e-4,
  :thickness => 0.0075)

laminateProperties =@compat Dict(
  :nplies => 5,
  :repeats => 1,
  :symmetric => false,
  :orientation => [0, 45, 90, 45, -45],
  :deltaTemp => -270.0,
  :deltaMoisture => 101
)

createLaminate!(laminateProperties, materialProperties, forcesAndMoments)
fm = [158.39, -41.6097, 210.166, 4.84751, 3.15249, 0.84751]
@assert(round(laminateProperties[:fTotal], 3) == round(fm, 3))
println("Tests for CLT passed!")
