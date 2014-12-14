abstract LaserElement

typealias LaserScheme Array{LaserElement, 1}

typealias FFTFun Union(Nothing, Function)

# JonesMatrix =================================================================
immutable type JonesMatrix{T<:Number} <: LaserElement
    m::Array{T, 2}
    
    function JonesMatrix(a::Array{T, 2})
        size(a) == (2,2) || error("invalid Jones matrix dimensions")
        new(a)
    end
end

# inner constructor requires explicit type declaration, outer one does not
JonesMatrix{T<:Number}(a::Array{T,2}) = JonesMatrix{T}(a)

*(M1::JonesMatrix, M2::JonesMatrix) = JonesMatrix(M1.m * M2.m)

function apply_Jones_matrix!(M::JonesMatrix, uX, uY)
    a = M.m
    for i = 1:length(uX)
        uX_ = uX[i]                 
        uY_ = uY[i]
        uX[i] = a[1,1] * uX_ + a[1,2] * uY_
        uY[i] = a[2,1] * uX_ + a[2,2] * uY_
    end
end


function Rotation(a=0)
    JonesMatrix([cos(a)  -sin(a); sin(a) cos(a)])
end

function Polarizer(a=0)
    polarizer = JonesMatrix([(1. + 0.im) 0; 0 0])
    Rotation(a) * polarizer * Rotation(-a)
end

function ArbitraryWavePlate(phase_shift, a=0.)
    dphi = phase_shift / 2.
    # signs are chosen so that X is fast axis and Y component
    # is slowed relatively to it
    plate = JonesMatrix([exp(-1.im*dphi) 0; 0 exp(1.im*dphi)])
    Rotation(a) * plate * Rotation(-a)
end

HalfWavePlate(a=0.) = ArbitraryWavePlate(pi, a)

QuarterWavePlate(a=0.) = ArbitraryWavePlate(pi/2, a)

# Fiber =======================================================================
immutable type Fiber{T<:Real} <:LaserElement
    L::T
    alpha::T
    betha::Vector{T} # maybe add separate bethas for both polarizations
    dbetha::T
    gamma::T
    gain::T
    gain_bandwidth::T
    saturation_energy::T
    max_steps::Integer
    adaptive_step::Bool
end

# no birefrigence case, dbetha = 0
Fiber{T<:Real}(L::T, alpha::T, betha::Vector{T}, gamma::T,
               gain::T, gain_bandwidth::T, saturation_energy::T,
               max_steps=1000::Integer, adaptive_step=false::Bool) =
    Fiber(L, alpha, betha, zero(T), gamma, gain, gain_bandwidth, saturation_energy,
          max_steps, adaptive_step)
# no gain case
FiberPassive{T<:Real}(L::T, alpha::T, betha::Vector{T}, dbetha::T, gamma::T,
               max_steps=1000::Integer, adaptive_step=false::Bool) =
    Fiber(L, alpha, betha, dbetha, gamma, zero(T), 1.e40*one(T), 1.e40*one(T),
          max_steps, adaptive_step)
# no gain and birefrigence case
FiberPassive{T<:Real}(L::T, alpha::T, betha::Vector{T}, gamma::T,
               max_steps=1000::Integer, adaptive_step=false::Bool) =
    FiberPassive(L, alpha, betha, zero(T), gamma, max_steps, adaptive_step)


# Pulse =======================================================================
type Pulse{Ty<:Real}
    uX::Vector{Complex{Ty}}
    uY::Vector{Complex{Ty}}
    t::Vector{Ty}
    w::Vector{Ty}
    fft_plan!::Function
    ifft_plan!::Function
end

apply_Jones_matrix!(M::JonesMatrix, p::Pulse) = apply_Jones_matrix!(M, p.uX, p.uY)

function +{T}(p1::Pulse{T}, p2::Pulse{T})
    (length(p1.t) == length(p2.t))              || error("unequal pulse lengths")
    # isequal(p1.fft_plan!, p1.fft_plan!)         || error("unequal FFT plans")
    # isequal(p1.ifft_plan!, p1.ifft_plan!)       || error("unequal IFFT plans")
    _T1 = calc_T(p1.t)
    _T2 = calc_T(p2.t)
    abs(_T1 - _T2) / middle(_T1, _T2) < 1.e-10   || error("unequal t grids")
    Pulse(p1.uX + p2.uX, p1.uY + p2.uY, p1.t, p1.w, p1.fft_plan!, p1.ifft_plan!)
end

function Pulse{T<:Real}(uX::Vector{Complex{T}}, uY::Vector{Complex{T}},
                        t::Vector{T}, w::Vector{T}, 
                        fft_plan! = nothing::FFTFun, ifft_plan! = nothing::FFTFun)
    (length(t) == length(w) == length(uX) == length(uY)) || 
        error("dimensions of all pulse arrays must match")

    u = similar(uX)
    fft_plan! == nothing && (fft_plan! = plan_fft!(u, (1,), FFTW.MEASURE))
    ifft_plan! == nothing && (ifft_plan! = plan_ifft!(u, (1,), FFTW.MEASURE))
    Pulse(copy(uX), copy(uY), copy(t), copy(w), fft_plan!, ifft_plan!)
end

function Pulse{T<:Real}(shape::Integer, T0::T, P0::T, C0::T, theta::T,
                        t::Vector{T}, t_offset::T, w::Vector{T},
                        fft_plan! = nothing::FFTFun, ifft_plan! = nothing::FFTFun)
    uX, uY = pulse_vec(shape, T0, P0, C0, theta, t, t_offset)
    Pulse(uX, uY, t, w, fft_plan!, ifft_plan!)
end

function Pulse{T<:Real}(shape::Integer, T0::T, P0::T, C0::T, theta::T,
                        n::Integer, T_::T, t_offset::T,
                        fft_plan! = nothing::FFTFun, ifft_plan! = nothing::FFTFun)
    t = t_grid(n, T_)
    w = w_grid(n, T_)
    uX, uY = pulse_vec(shape, T0, P0, C0, theta, t, t_offset)
    Pulse(uX, uY, t, w, fft_plan!, ifft_plan!)
end

# zero default time offset
Pulse{T<:Real}(shape::Integer, T0::T, P0::T, C0::T, theta::T, 
               t::Vector{T}, w::Vector{T},
               fft_plan! = nothing::FFTFun, ifft_plan! = nothing::FFTFun) =
    Pulse(shape, T0, P0, C0, theta, t, zero(T), w, fft_plan!, ifft_plan!)

Pulse{T<:Real}(shape::Integer, T0::T, P0::T, C0::T, theta::T,
               n::Integer, T_::T,
               fft_plan! = nothing::FFTFun, ifft_plan! = nothing::FFTFun) =
    Pulse(shape, T0, P0, C0, theta, n, T_, zero(T), fft_plan!, ifft_plan!)

function NoisePulse{T<:Real}(power::T, scale::T, t::Vector{T}, w::Vector{T},
                             fft_plan! = nothing::FFTFun, ifft_plan! = nothing::FFTFun)
    n = length(t)
    uX = zeros(Complex{T}, n);                  uY = zeros(Complex{T}, n)

    # create interpolation object
    T_ = calc_T(t)
    interp_t_grid = -T_:scale:T_
    ni = length(interp_t_grid)
    interp_u_grid = zeros(Complex{T}, ni)
    interp_theta_grid = zeros(T, ni)

    for i = 1:ni
        uabs = power * rand()
        phi = 2pi* rand()
        interp_u_grid[i]  = exp(1.im*phi) * uabs
        interp_theta_grid[i] = 2pi * rand()
    end

    interp_u_real = Grid.CoordInterpGrid(interp_t_grid, real(interp_u_grid),
                                     Grid.BCperiodic, Grid.InterpQuadratic)
    interp_u_imag = Grid.CoordInterpGrid(interp_t_grid, imag(interp_u_grid),
                                     Grid.BCperiodic, Grid.InterpQuadratic)
    interp_theta = Grid.CoordInterpGrid(interp_t_grid, interp_theta_grid,
                                     Grid.BCperiodic, Grid.InterpQuadratic)

    # interpolate pulse between random seed values
    for i = 1:n
        tc = t[i]
        u = interp_u_real[tc] + 1.im*interp_u_imag[tc]
        theta = interp_theta[tc]
        uX[i] = u * cos(theta)
        uY[i] = u * sin(theta)
    end

    Pulse(uX, uY, t, w, fft_plan!, ifft_plan!)
end

function NoisePulse{T<:Real}(power::T, scale::T, n::Integer, T_::T,
                             fft_plan! = nothing::FFTFun, ifft_plan! = nothing::FFTFun)
    t = t_grid(n, T_)
    w = w_grid(n, T_)
    NoisePulse(power, scale, t, w, fft_plan!, ifft_plan!)
end

function WhiteNoisePulse{T<:Real}(power::T, t::Vector{T}, w::Vector{T},
                                  fft_plan! = nothing::FFTFun, ifft_plan! = nothing::FFTFun)
    n = length(t)
    uX = zeros(Complex{T}, n);          uY = zeros(Complex{T}, n)

    u = similar(uX)
    fft_plan! == nothing && (fft_plan! = plan_fft!(u, (1,), FFTW.MEASURE))
    ifft_plan! == nothing && (ifft_plan! = plan_ifft!(u, (1,), FFTW.MEASURE))

    for i = 1:n
        uabs = power * rand()
        phi = 2pi* rand()
        u  = exp(1.im*phi) * uabs
        theta = 2pi * rand()

        uX[i] = u * cos(theta);         uY[i] = u * sin(theta)
    end
    fft_plan!(uX);                      fft_plan!(uY)              

    Pulse(uX, uY, t, w, fft_plan!, ifft_plan!)
end

function WhiteNoisePulse{T<:Real}(power::T, n::Integer, T_::T,
                                  fft_plan! = nothing::FFTFun, ifft_plan! = nothing::FFTFun)
    t = t_grid(n, T_)
    w = w_grid(n, T_)
    WhiteNoisePulse(power, t, w, fft_plan!, ifft_plan!)
end

# FileOutput ==================================================================
ONLY_X = true

type FileOutput <:LaserElement
    outdir::String
    postfix::String
    iteration::Integer
    onlyX::Bool
    first_iteration::Bool
    UX       # temporaries array
    UY       # temporaries array
end

FileOutput(outdir::String, onlyX=false::Bool) = 
    FileOutput(outdir, "", 1, onlyX, true, [0.im], [0.im])

FileOutput(outdir::String, postfix::String, onlyX=false::Bool) = 
    FileOutput(outdir, postfix, 1, onlyX, true, [0.im], [0.im])

# other elements ==============================================================
type RectangularSpectralFilter <: LaserElement
    bandwidth_ratio::Real
end

RectangularSpectralFilter() = RectangularSpectralFilter(0.8)

type GaussianSpectralFilter <: LaserElement
    bandwidth_fr::Real
end

GaussianSpectralFilter(wl0, bandwidth_wl) = 
    GaussianSpectralFilter(bandwidth_wl2fr_derivative(wl0, bandwidth_wl))

type SaturableAbsorber <: LaserElement
    modulation_depth::Real
    saturation_power::Real
    unsaturable_loss::Real
end

SaturableAbsorber(md, sp) = SaturableAbsorber(md, sp, 0.)

type SESAM <: LaserElement
    modulation_depth::Real
    saturation_power::Real
    recovery_time::Real
    unsaturable_loss::Real
end

SESAM(md, se, rt) = SESAM(md, se, rt, 0.)

type Coupler <: LaserElement
    transmittance
end

# PulseSensor =================================================================
immutable type PulseSensor <: LaserElement
    name::String
    # reported pulse parameters list    
end

PulseSensor() = PulseSensor("")

# ConvergenceDetector =========================================================
# TODO ...


