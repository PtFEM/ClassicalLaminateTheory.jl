import Base: show, showcompact

# Implementation of Classical Laminate Theory

function qmat(p::Dict)
  ν21 = p[:ν12]*p[:E2]/p[:E1]
  [p[:E1]/(1 - p[:ν12]*ν21) p[:ν12]*p[:E2]/(1 - p[:ν12]*ν21) 0;
  p[:ν12]*p[:E2]/(1 - p[:ν12]*ν21) p[:E2]/(1 - p[:ν12]*ν21) 0;
  0 0 p[:G12]]
end

function qbarmat(qmat::Array{Float64, 2}, theta::Float64)
  #theta = 2*theta*pi/360.0
  ct = cos(theta); st = sin(theta)
  c = cumprod([ct for i in 1:4], 1)
  s = cumprod([st for i in 1:4], 1)
  qbar = zeros(Float64, 3, 3)
  qbar[1, 1] = qmat[1, 1]*cos(theta)^4 + 2*(qmat[1, 2] + 2qmat[3, 3])*sin(theta)^2*cos(theta)^2 +
    qmat[2, 2]*sin(theta)^4
  qbar[1, 2] = (qmat[1, 1] + qmat[2, 2] - 4qmat[3, 3])*sin(theta)^2*cos(theta)^2 +
    qmat[1, 2]*(sin(theta)^4 + cos(theta)^4)
  qbar[2, 1] = qbar[1, 2]
  qbar[2, 2] = qmat[1, 1]*sin(theta)^4 + 2*(qmat[1, 2] + 2qmat[3, 3])*sin(theta)^2*cos(theta)^2 +
    qmat[2, 2]*cos(theta)^4
  qbar[1, 3] = (qmat[1, 1] - qmat[1, 2] - 2qmat[3, 3])*sin(theta)*cos(theta)^3 +
    (qmat[1, 2] - qmat[2, 2] + 2qmat[3, 3])*sin(theta)^3*cos(theta)
  qbar[3, 1] = qbar[1, 3]
  qbar[2, 3] = (qmat[1, 1] - qmat[1, 2] - 2qmat[3, 3])*sin(theta)^3*cos(theta) +
    (qmat[1, 2] - qmat[2, 2] + 2qmat[3, 3])*sin(theta)*cos(theta)^3
  qbar[3, 2] = qbar[2, 3]
  qbar[3, 3] = (qmat[1, 1] + qmat[2, 2] - 2qmat[1, 2] -2qmat[3, 3])*sin(theta)^2*cos(theta)^2 +
    qmat[3, 3]*(sin(theta)^4 + cos(theta)^4)
  qbar
end

function createLaminate!(l::Dict, m::Dict, f::Dict)
  l[:nLamina] = l[:nplies]*l[:repeats]*(l[:symmetric] ? 2 : 1)
  l[:laminateThickness] = l[:nLamina] * l[:tk]
  l[:repeatOrientation] = repeat(l[:orientation], outer=[l[:repeats]])
  l[:layerOrientation] = l[:symmetric] ? vcat(l[:repeatOrientation], reverse(l[:repeatOrientation])) : l[:orientation]
  if l[:symmetric]
    l[:z] = linspace(-l[:laminateThickness]/2, 0, l[:nplies]*l[:repeats]+1)
    l[:z] = l[:symmetric] ? vcat(l[:z], abs(reverse(l[:z]))[2:end]) : l[:z]
  else
    l[:z] = linspace(-l[:laminateThickness]/2, l[:laminateThickness]/2, l[:nLamina]+1)
  end
  computeABD!(l, m)
  l[:fTotal] = [f[:Nx],f[:Ny],f[:Nxy],f[:Mx],f[:My],f[:Mxy]] + l[:fT] + l[:fM]
  l[:f] = round(l[:ABD]\l[:fTotal], 10)
  l
end

function createLaminate!(l::Dict, m::Dict)
  l[:nLamina] = l[:nplies]*l[:repeats]*(l[:symmetric] ? 2 : 1)
  l[:laminateThickness] = l[:nLamina] * l[:tk]
  l[:repeatOrientation] = repeat(l[:orientation], outer=[l[:repeats]])
  l[:layerOrientation] = l[:symmetric] ? vcat(l[:repeatOrientation], reverse(l[:repeatOrientation])) : l[:orientation]
  if l[:symmetric]
    l[:z] = linspace(-l[:laminateThickness]/2, 0, l[:nplies]*l[:repeats]+1)
    l[:z] = l[:symmetric] ? vcat(l[:z], abs(reverse(l[:z]))[2:end]) : l[:z]
  else
    l[:z] = linspace(-l[:laminateThickness]/2, l[:laminateThickness]/2, l[:nLamina]+1)
  end
  
  computeABD!(l, m)
  l[:fTotal] = l[:fT] + l[:fM]
  l[:f] = round(l[:ABD]\l[:fTotal], 10)
  l
end

function computeABD!(l::Dict, m::Dict)
  A = zeros(Float64, 3, 3)
  B = zeros(Float64, 3, 3)
  D = zeros(Float64, 3, 3)
  ntemp = zeros(Float64, 3)
  mtemp = zeros(Float64, 3)
  nmoist = zeros(Float64, 3)
  mmoist = zeros(Float64, 3)
  for i in 1:l[:nLamina]
    tk = l[:z][i+1]-l[:z][i]
    tk2 = l[:z][i+1]^2-l[:z][i]^2
    tk3 = l[:z][i+1]^3-l[:z][i]^3
    
    theta = float(l[:layerOrientation][i])*(2*pi/360.0)
    A += qbarmat(qmat(m), theta)*tk
    B += (1/2) * qbarmat(qmat(m), theta)*tk2
    D += (1/3) * qbarmat(qmat(m), theta)*tk3
    
    alpha = [
      m[:alpha1]*cos(theta)^2 + m[:alpha2]*sin(theta)^2,
      m[:alpha1]*sin(theta)^2 + m[:alpha2]*cos(theta)^2,
      2*cos(theta)*sin(theta)*(m[:alpha1]-m[:alpha2])
    ]
    
    beta = [
      m[:beta1]*cos(theta)^2 + m[:beta2]*sin(theta)^2,
      m[:beta1]*sin(theta)^2 + m[:beta2]*cos(theta)^2,
      2*cos(theta)*sin(theta)*(m[:beta1]-m[:beta2])
    ]
    
    ntemp += qbarmat(qmat(m), theta)*alpha*l[:deltaTemp]*tk
    mtemp += (1/2)*qbarmat(qmat(m), theta)*beta/100.0*l[:deltaTemp]*tk2
    nmoist += qbarmat(qmat(m), theta)*alpha*l[:deltaMoisture]*tk
    mmoist += (1/2)*qbarmat(qmat(m), theta)*beta/100.0*l[:deltaMoisture]*tk2
  end
  
  l[:fT] = vcat(ntemp, mtemp)
  l[:fM] = vcat(nmoist, mmoist)
  l[:ABD] = [ A B; B D]
  l[:abd] = inv(l[:ABD])
  l[:Ex] = 1/(l[:laminateThickness]*l[:abd][1, 1])
  l[:Ey] = 1/(l[:laminateThickness]*l[:abd][2, 2])
  l[:Gxy] = 1/(l[:laminateThickness]*l[:abd][3, 3])
  l[:νxy] = -l[:abd][2, 1] / l[:abd][1, 1]
  l[:Eᶠx] = 12/(l[:laminateThickness]^3*l[:abd][4, 4])
  l[:Eᶠy] = 12/(l[:laminateThickness]^3*l[:abd][5, 5])
  l[:Gᶠxy] = 12/(l[:laminateThickness]^3*l[:abd][6, 6])
  l[:νᶠxy] = -l[:abd][5, 4] / l[:abd][4, 4]
  l[:νᶠyx] = -l[:abd][5, 4] / l[:abd][5, 5]
end


function model_show(io::IO, m::Dict, compact::Bool)
  if compact
    println(m)
  else
    println("Effective material Properties:")
    println("  Ex  (psi) =       $(m[:Ex])")
    println("  Ey  (psi) =       $(m[:Ey])")
    println("  Gxy (psi) =       $(m[:Gxy])")
    println("  νxy (psi) =       $(m[:νxy])")
    println()
    println("Strain Curvature (EpsilonKappa):")
    println("  εx   =       $(m[:f][1])")
    println("  εy   =       $(m[:f][2])")
    println("  γxy  =       $(m[:f][3])")
    println("  κx   =       $(m[:f][4])")
    println("  κy   =       $(m[:f][5])")
    println("  κxy  =       $(m[:f][6])")
    println()
    println("Effective flexural Properties:")
    println("  Eᶠx  (psi) =       $(m[:Eᶠx])")
    println("  Eᶠy  (psi) =       $(m[:Eᶠy])")
    println("  Gᶠxy (psi) =       $(m[:Gᶠxy])")
    println("  νᶠxy (psi) =       $(m[:νᶠxy])")
    println("  νᶠyx (psi) =       $(m[:νᶠyx])")
    println()
  end
end

show(io::IO, m::Dict) = model_show(io, m, false)
showcompact(io::IO, m::Dict) = model_show(io, m, true)

