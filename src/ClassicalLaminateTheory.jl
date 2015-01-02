module ClassicalLaminateTheory

# package code goes here
using Compat

include("clt.jl")

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
  
  # From ClassicalLaminateTheory.jl
  JULIA_SVG_BROWSER

end # module
