abstract LaserElement

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
Fiber(alpha, betha, dbetha, gamma) = Fiber(alpha, betha, dbetha, gamma, 0., 1.e40, 1.e40)

# no birefrigence case
Fiber(alpha, betha, gamma) = Fiber(alpha, betha, 0., gamma)

type Pulse{T<:Real}
    uX::Vector{Complex{T}}
    uY::Vector{Complex{T}}
    t::Vector{T}
    w::Vector{T}
    fft_plan!::Function
    ifft_plan!::Function
end

function Pulse(uX, uY, t, w)
    u = copy(uX)
    fft_plan! = plan_fft!(u, (1,), FFTW.MEASURE)
    ifft_plan! = plan_ifft!(u, (1,), FFTW.MEASURE)
    Pulse(uX, uY, t, w, fft_plan!, ifft_plan!)
end    

type FileOutput <:LaserElement
    outdir::String
    iteration::Integer
end

FileOutput(outdir::String) = FileOutput(outdir, 0)

propagate_through!(p::Pulse, M::JonesMatrix) = apply_Jones_matrix(M, p.uX, p.uY)

function propagate_through!(p::Pulse, fiber::Fiber)
    u_plotX, u_plotY, U_plotX, U_plotY, n_steps, n_steps_rejected, steps =
        rk4ip_vec!(p, fiber)    
end

function propagate_through!(p::Pulse, fout::FileOutput)
    UX = similar(p.uX)
    UY = similar(p.uY)
    T = (p.t[end] - p.t[1]) / 2   
    spectrum!(p.uX, UX, p.ifft_plan!, T)
    spectrum!(p.uY, UY, p.ifft_plan!, T)

    i = fout.iteration
    fwrite(outdir, "uX", i, p.uX)
    fwrite(outdir, "uY", i, p.uY)
    fwrite(outdir, "UX", i, UX)
    fwrite(outdir, "UY", i, UY)

    fout.iteration += 1
end

function run_laser(n, T_window, scheme, L, T0, P0, C0, theta, shape=0)
    P1 = W4(a1)
    P2 = W4(a2)
    P3 = W2(a3)
    P23p = P2 * P3 * [1 0; 0 0]

    uX, uY = pulse_vec(shape, T0, P0, C0, theta, t)

    outdir = mkpath_today()
    
    n_iter = 10

    for i in 1:n_iter
        apply_Jones_matrix!(P1, uX, uY)

        u_plotX, u_plotY, U_plotX, U_plotY, n_steps, n_steps_rejected, steps =
            rk4ip_vec!(uX, uY, L, 1.e-10L, t, w, alpha, betha, dbetha, gamma,
                       gain, gain_bandwidth, E_sat, fft_plan!, ifft_plan!, 0, 0)

        apply_Jones_matrix!(P23p, uX, uY)

        spectrum!(uX, UX, ifft_plan!, T)
        spectrum!(uY, UY, ifft_plan!, T)

        fwrite(outdir, "uX", i, uX)
        fwrite(outdir, "uY", i, uY)
        fwrite(outdir, "UX", i, UX)
        fwrite(outdir, "UY", i, UY)       
    end
end

