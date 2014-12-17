module NLSE

using NumericExtensions, NumericFuns, Devectorize, Dates, Grid

include("laser_scheme.jl")
include("laser.jl")

include("output.jl")
include("utils.jl")
include("run.jl")
include("rk4ip_utils.jl")

include("rk4ip.jl")
include("test_cases.jl")

include("rk4ip_vec.jl")
include("test_cases_vec.jl")

include("ssfm.jl")
end


