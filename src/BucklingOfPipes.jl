module BucklingOfPipes

using Compat

include("ClassicalLaminateTheory:,"clt.jl")
include("basisfunctions.jl")

if !isdefined(Main, :JULIA_SVG_BROWSER)
  JULIA_SVG_BROWSER = ""
  try
    JULIA_SVG_BROWSER = ENV["JULIA_SVG_BROWSER"]
  catch e
    println("Environment variable JULIA_SVG_BROWSER not found.")
    JULIA_SVG_BROWSER = ""
  end
end

export
  # From clt.jl
  qmat,
  qbarmat,
  createLaminate!,
  
  # From basisfunctions
  basisfunctions,
  basisfunction,
  
  # From BucklingOfPipes.jl
  JULIA_SVG_BROWSER

end # module
