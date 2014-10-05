module NLSE

using Lumberjack, NumericExtensions, NumericFuns, Devectorize

export t_grid, w_grid, secant_pulse, gaussian_pulse, integrate_RK4IP
export x, y

include("run.jl")
include("utils.jl")
include("ssfm.jl")
include("rk4ip.jl")

end


