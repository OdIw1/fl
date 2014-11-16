abstract LaserElement

typealias LaserScheme Array{LaserElement, 1}

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

function Rotation(a=0)
    JonesMatrix([cos(a)  -sin(a); sin(a) -cos(a)])
end

function Polarizer(a=0)
    polarizer = JonesMatrix([(1. + 0.im) 0; 0 0])
    Rotation(a) * polarizer * Rotation(-a)
end

function HalfWavePlate(a=0)
    half_wave = JonesMatrix([-1.im 0; 0 1.im])
    Rotation(a) * half_wave * Rotation(-a)
end

function QuarterWavePlate(a=0.)
    quarter_wave = JonesMatrix(1. / sqrt(2) * [(1. - 1.im)  0; 0 (1. + 1.im)])
    Rotation(a) * quarter_wave * Rotation(-a)
end

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
end

# no birefrigence case
Fiber{T}(L::T, alpha::T, betha::Vector{T}, gamma::T,
         gain::T, gain_bandwidth::T, saturation_energy::T) =
    Fiber(L, alpha, betha, 0., gamma, gain, gain_bandwidth, saturation_energy)

# no gain case
Fiber{T}(L::T, alpha::T, betha::Vector{T}, dbetha::T, gamma::T) =
    Fiber(L, alpha, betha, dbetha, gamma, 0., 1.e40, 1.e40)
# no gain and birefrigence case
Fiber{T}(L::T, alpha::T, betha::Vector{T}, gamma::T) = Fiber(L, alpha, betha, 0., gamma)

# FileOutput ==================================================================
type FileOutput <:LaserElement
    outdir::String
    postfix::String
    iteration::Integer
end

FileOutput(outdir::String) = FileOutput(outdir, "", 1)

FileOutput(outdir::String, postfix::String) = FileOutput(outdir, postfix, 1)

# Pulse =======================================================================
type Pulse{Ty<:Real}
    uX::Vector{Complex{Ty}}
    uY::Vector{Complex{Ty}}
    t::Vector{Ty}
    w::Vector{Ty}
    n::Integer
    T::Ty
    fft_plan!::Function
    ifft_plan!::Function
end

function Pulse{T}(uX::Vector{Complex{T}}, uY::Vector{Complex{T}},
                  t::Vector{T}, w::Vector{T}, 
                  fft_plan! = nothing::Union(Nothing, Function),
                  ifft_plan! = nothing::Union(Nothing, Function))
    n = length(t)
    dt = (t[end] - t[1]) / (n - 1)   
    _T = (dt + t[end] - t[1]) / 2

    (n == length(w) == length(uX) == length(uY)) || error ("dimensions of all arrays must match")

    u = copy(uX)
    fft_plan! == nothing && (fft_plan! = plan_fft!(u, (1,), FFTW.MEASURE))
    ifft_plan! == nothing && (ifft_plan! = plan_ifft!(u, (1,), FFTW.MEASURE))
    Pulse(uX, uY, t, w, n, _T, fft_plan!, ifft_plan!)
end

function Pulse{T}(shape::Integer, T0::T, P0::T, C0::T, theta::T,
                  t::Vector{T}, t_offset::T, w::Vector{T},
                  fft_plan! = nothing::Union(Nothing, Function),
                  ifft_plan! = nothing::Union(Nothing, Function))
    uX, uY = pulse_vec(shape, T0, P0, C0, theta, t, t_offset)
    Pulse(uX, uY, t, w, fft_plan!, ifft_plan!)
end

# zero default time offset
Pulse{T}(shape::Integer, T0::T, P0::T, C0::T, theta::T, t::Vector{T}, w::Vector{T},
         fft_plan! = nothing::Union(Nothing, Function),
         ifft_plan! = nothing::Union(Nothing, Function)) =
    Pulse(shape, T0, P0, C0, theta, t, 0., w, fft_plan!, ifft_plan!)

# secant pulse with zero offset
Pulse{T}(T0::T, P0::T, C0::T, theta::T, t::Vector{T}, w::Vector{T},
         fft_plan! = nothing::Union(Nothing, Function),
         ifft_plan! = nothing::Union(Nothing, Function)) =
    Pulse(0, T0, P0, C0, theta, t, 0., w, fft_plan!, ifft_plan!)

# ConvergenceDetector =========================================================
# TODO ...


function apply_Jones_matrix!(M::JonesMatrix, uX, uY)
    a = M.m
    for i = 1:length(uX)
        uX_ = uX[i]                 
        uY_ = uY[i]
        uX[i] = a[1,1] * uX_ + a[1,2] * uY_
        uY[i] = a[2,1] * uX_ + a[2,2] * uY_
    end
end

apply_Jones_matrix!(M::JonesMatrix, p::Pulse) = apply_Jones_matrix!(M, p.uX, p.uY)
