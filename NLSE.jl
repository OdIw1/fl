module NLSE

using NumericExtensions, NumericFuns, Devectorize, Datetime

include("output.jl")
include("utils.jl")
include("run.jl")

include("rk4ip.jl")
include("test_cases.jl")

include("rk4ip_vec.jl")
include("test_cases_vec.jl")

include("polarization_control.jl")

include("ssfm.jl")

end


