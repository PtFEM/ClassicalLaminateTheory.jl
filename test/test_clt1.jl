using Compat, ClassicalLaminateTheory

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
  :beta2 => 20e-4,
  :thickness => 0.0075)

laminateProperties = @compat Dict(
  :nplies => 4,
  :repeats => 1,
  :symmetric => true,
  :orientation => [0, 45, 90, -45],
  :deltaTemp => 0.0,
  :deltaMoisture => 0.00
)

createLaminate!(laminateProperties, materialProperties, forcesAndMoments)

ABDmat = [
  525698.0  155811.0       0.0    0.0    -0.0     0.0;
  155811.0  525698.0       0.0   -0.0     0.0     0.0;
       0.0       0.0  184943.0    0.0     0.0     0.0;
       0.0      -0.0       0.0  250.495  32.687  23.619;
      -0.0       0.0       0.0   32.687  93.036  23.619;
       0.0       0.0       0.0   23.619  23.619  41.427
];

@assert round(laminateProperties[:ABD]) == round(ABDmat)
println("Tests for CLT1 passed!")
