abstract LaserElement

typealias LaserScheme Array{LaserElement, 1}

immutable type JonesMatrix{T<:Number} <: LaserElement
    m::Array{T, 2}

    function JonesMatrix(a)
        isa(a, Array{T, 2}) || error("Jones matrix must be an array")
        size(a) == (2,2) || error("invalid Jones matrix dimensions")
        new(m)
    end
end

*(M1::JonesMatrix, M2::JonesMatrix) = JonesMatrix(M1.m * M2.m)

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

# no gain case
Fiber{T}(L::T, alpha::T, betha::Vector{T}, dbetha::T, gamma::T) = 
    Fiber(L, alpha, betha, dbetha, gamma, 0., 1.e40, 1.e40)
# no birefrigence case
Fiber{T}(L::T, alpha::T, betha::Vector{T}, gamma::T) = Fiber(L, alpha, betha, 0., gamma)

type FileOutput <:LaserElement
    outdir::String
    postfix::String
    iteration::Integer
end

FileOutput(outdir::String) = FileOutput(outdir, "")

FileOutput(outdir::String, postfix::String) = FileOutput(outdir, postfix, 0)

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
                  t::Vector{T}, w::Vector{T})
    n = length(t)
    dt = (t[end] - t[1]) / (n - 1)   
    T = (dt + t[end] - t[1]) / 2

    (n == length(w) == length(uX) == length(uY)) || error ("dimensions of all arrays must match")

    u = copy(uX)
    fft_plan! = plan_fft!(u, (1,), FFTW.MEASURE)
    ifft_plan! = plan_ifft!(u, (1,), FFTW.MEASURE)
    Pulse(uX, uY, t, w, n, T, fft_plan!, ifft_plan!)
end

function Pulse{T}(shape::Integer, T0::T, P0::T, C0::T,
                  t::Vector{T}, t_offset::T, w::Vector{T})
    uX, uY = pulse_vec(shape, T0, P0, C0, t, t_offset)
    Pulse(uX, uY, t, w)
end

# similar, but with zero default time offset
function Pulse{T}(shape::Integer, T0::T, P0::T, C0::T,
                  t::Vector{T}, w::Vector{T})
    uX, uY = pulse_vec(shape, T0, P0, C0, t, 0.)
    Pulse(uX, uY, t, w)
end


    

