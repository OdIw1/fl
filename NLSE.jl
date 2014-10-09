module NLSE

using Lumberjack, NumericExtensions, NumericFuns, Devectorize

include("rk4ip.jl")
include("ssfm.jl")
include("utils.jl")
include("output.jl")
include("run.jl")
include("test_cases.jl")

end


