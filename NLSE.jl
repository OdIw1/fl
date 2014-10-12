module NLSE

using NumericExtensions, NumericFuns, Devectorize

include("rk4ip.jl")
include("test_cases.jl")

include("rk4ip_vec.jl")
include("test_cases_vec.jl")

include("ssfm.jl")
include("utils.jl")
include("output.jl")
include("run.jl")

end


