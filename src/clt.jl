import Base: show, showcompact

# Implementation of Classical Laminate Theory

function qmat(m::Dict, mat::Int64=1)
  res = zeros(3, 3)
  if typeof(m[:E1]) == Float64
    @assert mat == 1
    m[:ν21] = m[:ν12]*m[:E2]/m[:E1]
    res = [m[:E1]/(1 - m[:ν12]*m[:ν21]) m[:ν12]*m[:E2]/(1 - m[:ν12]*m[:ν21]) 0;
    m[:ν12]*m[:E2]/(1 - m[:ν12]*m[:ν21]) m[:E2]/(1 - m[:ν12]*m[:ν21]) 0;
    0 0 m[:G12]]
  elseif typeof(m[:E1]) == Array{Float64, 1}
    m[:ν21] = similar(m[:ν12], Float64)
    for i in 1:length(m[:E1])
      m[:ν21][i] = m[:ν12][i]*m[:E2][i]/m[:E1][i]
    end
    res = 
      [m[:E1][mat]/(1 - m[:ν12][mat]*m[:ν21][mat]) m[:ν12][mat]*m[:E2][mat]/(1 - m[:ν12][mat]*m[:ν21][mat]) 0;
      m[:ν12][mat]*m[:E2][mat]/(1 - m[:ν12][mat]*m[:ν21][mat]) m[:E2][mat]/(1 - m[:ν12][mat]*m[:ν21][mat]) 0;
      0 0 m[:G12][mat]
    ]
  end
  res
end

function qbarmat(qmat::Array{Float64, 2}, theta::Float64)
  c1 = cos(theta)
  s1 = sin(theta)
  c2 = c1*c1
  s2 = s1*s1
  c3 = c1*c2
  s3 = s1*s2
  c4 = c2*c2
  s4 = s2*s2
  qbar = zeros(Float64, 3, 3)
  qbar[1, 1] = qmat[1, 1]*c4 + 2*(qmat[1, 2] + 2qmat[3, 3])*s2*c2 + qmat[2, 2]*s4
  qbar[1, 2] = (qmat[1, 1] + qmat[2, 2] - 4qmat[3, 3])*s2*c2 + qmat[1, 2]*(s4 + c4)
  qbar[2, 1] = qbar[1, 2]
  qbar[2, 2] = qmat[1, 1]*s4 + 2*(qmat[1, 2] + 2qmat[3, 3])*s2*c2 + qmat[2, 2]*c4
  qbar[1, 3] = (qmat[1, 1] - qmat[1, 2] - 2qmat[3, 3])*s1*c3 +
    (qmat[1, 2] - qmat[2, 2] + 2qmat[3, 3])*s3*c1
  qbar[3, 1] = qbar[1, 3]
  qbar[2, 3] = (qmat[1, 1] - qmat[1, 2] - 2qmat[3, 3])*s3*c1 +
    (qmat[1, 2] - qmat[2, 2] + 2qmat[3, 3])*s1*c3
  qbar[3, 2] = qbar[2, 3]
  qbar[3, 3] = (qmat[1, 1] + qmat[2, 2] - 2qmat[1, 2] - 2qmat[3, 3])*s2*c2 +
    qmat[3, 3]*(s4 + c4)
  qbar
end

function createLaminate!(l::Dict, m::Dict, f::Dict)
  l[:thickness] = 0.0
  if !(:materials in keys(l)) || (:materials in keys(l) && typeof(l[:materials]) == Int64)
    l[:materials] = ones(Int64, l[:nplies])
  end
  if :repeats in keys(l) && l[:repeats] > 1
    l[:nplies] *= l[:repeats]
    l[:orientation] = repeat(l[:orientation], outer=[l[:repeats]])
    l[:materials] = repeat(l[:materials], outer=[l[:repeats]])
  end
  if :symmetric in keys(l) && l[:symmetric] == true
    l[:nplies] *= 2
    l[:orientation] = vcat(l[:orientation], reverse(l[:orientation]))
    l[:materials] = vcat(l[:materials], reverse(l[:materials]))
  end
  if :symmetricrepeats in keys(l) && l[:symmetricrepeats] > 1
    l[:nplies] *= l[:symmetricrepeats]
    l[:orientation] = repeat(l[:orientation], outer=[l[:symmetricrepeats]])
    l[:materials] = repeat(l[:materials], outer=[l[:symmetricrepeats]])
  end
  l[:z] = zeros(Float64, l[:nplies]+1)
  if typeof(m[:thickness]) == Float64
    l[:thickness] = l[:nplies] * m[:thickness]
    l[:z] = linspace(-l[:thickness]/2, l[:thickness]/2, l[:nplies]+1)
  elseif typeof(m[:thickness]) == Array{Float64, 1}
    @assert l[:nplies] == length(l[:materials])
    for i in 1:l[:nplies]
      l[:thickness] += m[:thickness][l[:materials][i]]
    end
    l[:z][1] = -l[:thickness]/2.0
    for i in 1:l[:nplies]
      l[:z][i+1] = l[:z][i] + m[:thickness][l[:materials][i]]
    end
  end
  computeABD!(l, m)
  l[:fTotal] = [f[:Nx],f[:Ny],f[:Nxy],f[:Mx],f[:My],f[:Mxy]] + l[:fT] + l[:fM]
  l[:f] = round(l[:ABD]\l[:fTotal], 15)
end

function computeABD!(l::Dict, m::Dict)
  A = zeros(Float64, 3, 3)
  B = zeros(Float64, 3, 3)
  D = zeros(Float64, 3, 3)
  ntemp = zeros(Float64, 3)
  mtemp = zeros(Float64, 3)
  nmoist = zeros(Float64, 3)
  mmoist = zeros(Float64, 3)
  for i in 1:l[:nplies]
    tk = l[:z][i+1]-l[:z][i]
    tk2 = l[:z][i+1]^2-l[:z][i]^2
    tk3 = l[:z][i+1]^3-l[:z][i]^3
    
    theta = float(l[:orientation][i])*(pi/180.0)
    A += qbarmat(qmat(m, l[:materials][i]), theta)*tk
    B += (1/2) * qbarmat(qmat(m, l[:materials][i]), theta)*tk2
    D += (1/3) * qbarmat(qmat(m, l[:materials][i]), theta)*tk3
    
    alpha = [
      m[:alpha1][l[:materials][i]]*cos(theta)^2 + m[:alpha2][l[:materials][i]]*sin(theta)^2,
      m[:alpha1][l[:materials][i]]*sin(theta)^2 + m[:alpha2][l[:materials][i]]*cos(theta)^2,
      2*cos(theta)*sin(theta)*(m[:alpha1][l[:materials][i]]-m[:alpha2][l[:materials][i]])
    ]
    
    beta = [
      m[:beta1][l[:materials][i]]*cos(theta)^2 + m[:beta2][l[:materials][i]]*sin(theta)^2,
      m[:beta1][l[:materials][i]]*sin(theta)^2 + m[:beta2][l[:materials][i]]*cos(theta)^2,
      2*cos(theta)*sin(theta)*(m[:beta1][l[:materials][i]]-m[:beta2][l[:materials][i]])
    ]
    
    ntemp += qbarmat(qmat(m, l[:materials][i]), theta)*alpha*l[:deltaTemp]*tk
    mtemp += (1/2)*qbarmat(qmat(m, l[:materials][i]), theta)*beta*l[:deltaTemp]*tk2
    nmoist += qbarmat(qmat(m, l[:materials][i]), theta)*alpha*l[:deltaMoisture]*tk
    mmoist += (1/2)*qbarmat(qmat(m, l[:materials][i]), theta)*beta*l[:deltaMoisture]*tk2
  end
  
  l[:fT] = vcat(ntemp, mtemp/100.0)
  l[:fM] = vcat(nmoist*100.0, mmoist)
  l[:ABD] = [ A B; B D]
  l[:abd] = inv(l[:ABD])
  l[:Ex] = 1/(l[:thickness]*l[:abd][1, 1])
  l[:Ey] = 1/(l[:thickness]*l[:abd][2, 2])
  l[:Gxy] = 1/(l[:thickness]*l[:abd][3, 3])
  l[:νxy] = -l[:abd][2, 1] / l[:abd][1, 1]
  l[:Eᶠx] = 12/(l[:thickness]^3*l[:abd][4, 4])
  l[:Eᶠy] = 12/(l[:thickness]^3*l[:abd][5, 5])
  l[:Gᶠxy] = 12/(l[:thickness]^3*l[:abd][6, 6])
  l[:νᶠxy] = -l[:abd][5, 4] / l[:abd][4, 4]
  l[:νᶠyx] = -l[:abd][5, 4] / l[:abd][5, 5]
end


function model_show(io::IO, l::Dict, compact::Bool)
  if compact
    println(m)
  else
    println("Effective material Properties:")
    println("  Ex  (psi) =       $(l[:Ex])")
    println("  Ey  (psi) =       $(l[:Ey])")
    println("  Gxy (psi) =       $(l[:Gxy])")
    println("  νxy       =       $(l[:νxy])")
    #println("  νyx       =       $(l[:νyx])")
    println()
    println("Strain Curvature (EpsilonKappa):")
    println("  εx   =       $(l[:f][1])")
    println("  εy   =       $(l[:f][2])")
    println("  γxy  =       $(l[:f][3])")
    println("  κx   =       $(l[:f][4])")
    println("  κy   =       $(l[:f][5])")
    println("  κxy  =       $(l[:f][6])")
    println()
    println("Effective flexural Properties:")
    println("  Eᶠx  (psi) =       $(l[:Eᶠx])")
    println("  Eᶠy  (psi) =       $(l[:Eᶠy])")
    println("  Gᶠxy (psi) =       $(l[:Gᶠxy])")
    println("  νᶠxy       =       $(l[:νᶠxy])")
    println("  νᶠyx       =       $(l[:νᶠyx])")
    println()
  end
end

show(io::IO, m::Dict) = model_show(io, m, false)
showcompact(io::IO, m::Dict) = model_show(io, m, true)

