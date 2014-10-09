function beta_pskm_to_sm(beta...)
    b = zeros(Float64, length(beta))
    for k = 1:length(beta)
        b[k] = beta[k] * (1e-12)^(k+1) * 1.e-3
    end
    b
end

function beta_psm_to_sm(beta...)
    b = zeros(Float64, length(beta))
    for k = 1:length(beta)
        b[k] = beta[k] * (1e-12)^(k+1)
    end
    b
end

function pulse_propagation_params(T0, P0, gamma, beta...)
    ld = zeros(Float64, length(beta))
    soliton_order = zeros(ld)

    ln = 1. / (gamma * P0)
    for k = 1:length(beta)
        ld[k] = T0^(k+1) / abs(beta[k])
        soliton_order[k] = sqrt(ld[k] / ln)
    end
    ln, ld, soliton_order
end

function t_grid(n, T)
# In MATLAB version 'T' was a scaling factor: time window was (-T0 T, T0 T),
# where T0 is pulse FWHM;  in Julia time window is (-T, T)
# also 'n' is no longer binary log of the number of grid points, it is the number itself
    dt = 2 * T / n
    dt * [-n/2 : (n/2 - 1)]
end

function w_grid(n, T)
    dw = pi / T
    dw * [0 : (n/2 -1), -n/2 : -1]
end

function secant_pulse(T0, P0, C0, t::Vector, t_offset=0)
    t_scaled = (t - t_offset) / T0
    sqrt(P0) * sech(t_scaled) .* exp(-0.5im * C0 * (t_scaled).^2)
end

function gaussian_pulse(m_order, T0, P0, C0, t::Vector, t_offset=0)
    t_scaled = (t - t_offset) / T0
    sqrt(P0) * exp(-0.5(1 + 1im * C0) * t_scaled .^ (2m_order))
end


function spectrum(u, ifft_plan!, T)
    U = similar(u)
    spectrum!(u, U, ifft_plan!, T)
    U
end

function spectrum!(u, U, ifft_plan!, T)
    n = length(u)
    BLAS.blascopy!(n, u, 1, U, 1)
    ifft_plan!(U)
    fftshift!(U)
    BLAS.scal!(n, T * sqrt(2/pi), U, 1)
end

function fftshift!(u)
    n = length(u)
    if mod(n, 2) != 0 
        throw(ArgumentError("fftshift! is defined only for even-length arrays"))
    end

    shift = div(n, 2)
    for i = 1:shift
        u[i], u[shift+i] = u[shift+i], u[i]
    end
end

function df!(u, du, dx)
    # du = d(u)/dx, u and du MUST be different arrays
    du[1] = (u[2] - u[end]) / (2dx)
    du[end] = (u[1] - u[end-1]) / (2dx)
    @simd for i in 2:(length(u)-1)
        @inbounds du[i] = (u[i+1] - u[i-1]) / (2dx)
    end
end

function integration_error_local{T1<:Complex, T2<:Real}(u1::Vector{T1}, u2::Vector{T1},
                                                        ue_::Vector{T2}, atol=1.e-6, rtol=1.e-6)
    @simd for i in 1:length(u1)
        @inbounds ue_[i] = abs(u1[i] - u2[i]) / (atol + rtol * max(abs(u1[i]), abs(u2[i])))
    end
    maximum(ue_)
    # maximum(abs((abs(u1 - u2) ./ error_scale)));
end

function integration_error_global{T<:Complex}(u1::Vector{T}, u2::Vector{T},
                                              ue_cplx_::Vector{T}, atol=1.e-5, rtol=1.e-5)
    n = length(u1)
    BLAS.blascopy!(n, u2, 1, ue_cplx_, 1)
    BLAS.axpy!(n, -1. + 0im, u1, 1, ue_cplx_, 1)
    BLAS.nrm2(n, ue_cplx_, 1) / (atol + rtol * BLAS.nrm2(n, u2, 1))
end