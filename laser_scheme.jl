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
Fiber{T<:Real}(L::T, alpha::T, betha::Vector{T}, dbetha::T, gamma::T,
               max_steps=1000::Integer, adaptive_step=false::Bool) =
    Fiber(L, alpha, betha, dbetha, gamma, zero(T), 1.e40*one(T), 1.e40*one(T),
          max_steps, adaptive_step)
# no gain and birefrigence case
Fiber{T<:Real}(L::T, alpha::T, betha::Vector{T}, gamma::T,
               max_steps=1000::Integer, adaptive_step=false::Bool) =
    Fiber(L, alpha, betha, zero(T), gamma, max_steps, adaptive_step)

# FileOutput ==================================================================
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

# other elements ==============================================================
type SpectralFilter <: LaserElement
    bandwidth::Real
end

SpectralFilter() = SpectralFilter(0.8)

type SaturableAbsorber <: LaserElement
    modulation_depth
    saturation_power
end

type Coupler <: LaserElement
    transmittance
end

# PulseSensor =================================================================
immutable type PulseSensor <: LaserElement
    name::String
    # reported pulse parameters list    
end

PulseSensor() = PulseSensor("unnamed_sensor")

# ConvergenceDetector =========================================================
# TODO ...


