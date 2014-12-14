ADAPTIVE_STEP = true
FIXED_STEP = false

function pulse_energy(u, dt)
    n = length(u)
    energy = sqr(BLAS.nrm2(n, u, 1)) * dt    
end 

function gain_spectral_factor(w, g, g_bandwidth)
    0.5g ./ (1 + w.^2 / g_bandwidth^2)
end

function gain_saturation_vec(uX, uY, dt, saturation_energy)
    n = length(uX)
    EX = sqr(BLAS.nrm2(n, uX, 1)) * dt;             EY = sqr(BLAS.nrm2(n, uY, 1)) * dt
    energy = EX + EY
    1 / (1. + energy / saturation_energy)
end

function gain_saturation_scal(u, dt, saturation_energy)
    n = length(u)
    energy = sqr(BLAS.nrm2(n, u, 1)) * dt
    1 / (1. + energy / saturation_energy)
end

function dispersion_exp_vec!(d_exp, d_no_gain, g_spec_factor, uX, uY, dt, saturation_energy)
    g_sat = gain_saturation_vec(uX, uY, dt, saturation_energy)
    @devec d_exp[:] = d_no_gain + g_sat .* g_spec_factor
end

function dispersion_exp_scal!(d_exp, d_no_gain, g_spec_factor, u, dt, saturation_energy)
    g_sat = gain_saturation_scal(u, dt, saturation_energy)
    @devec d_exp[:] = d_no_gain + g_sat .* g_spec_factor
end

function D_exp_no_gain(w, alpha, betha)
    # pay attention to the order of Fourier transforms that determine
    # the sign of differentiation operator
    res = zeros(Complex{Float64}, length(w))
    res -= alpha/2
    for k in 1:length(betha)
        res += (1.im / factorial(k+1) * betha[k]) * w.^(k+1) # * (-1)^(k+1)
    end
    return res
end

function df!(u, du, dx)
    # du = d(u)/dx, u and du MUST be different arrays
    du[1] = (u[2] - u[end]) / (2dx)
    du[end] = (u[1] - u[end-1]) / (2dx)
    @simd for i in 2:(length(u)-1)
        @inbounds du[i] = (u[i+1] - u[i-1]) / (2dx)
    end
end

# atol should be adjusting according to pulse energy
function integration_error_local{T1<:Complex, T2<:Real}(u1::Vector{T1}, u2::Vector{T1},
                                                        ue_::Vector{T2}, atol=1.e-20, rtol=1.e-6)
    @simd for i in 1:length(u1)
        @inbounds ue_[i] = abs(u1[i] - u2[i]) / (atol + rtol * max(abs(u1[i]), abs(u2[i])))
    end
    maximum(ue_)
    # maximum(abs((abs(u1 - u2) ./ error_scale)));
end

function integration_error_global{T<:Complex}(u1::Vector{T}, u2::Vector{T},
                                              ue_cplx_::Vector{T}, atol=1.e-20, rtol=1.e-6)
    n = length(u1)
    BLAS.blascopy!(n, u2, 1, ue_cplx_, 1)
    BLAS.axpy!(n, -1. + 0im, u1, 1, ue_cplx_, 1)
    BLAS.nrm2(n, ue_cplx_, 1) / (atol + rtol * BLAS.nrm2(n, u2, 1))
end

function PI_control_factor(err, err_prev, ae=0.7, be=0.4)
    err^(-ae/5) * err_prev^(be/5)
end

function scale_step_fail(err, err_prev, ae=0.7, be=0.4)
    0.8max(1/5., PI_control_factor(err, err_prev, ae, be))
end

function scale_step_ok(err, err_prev, ae=0.7, be=0.4)
    0.8min(10, PI_control_factor(err, err_prev, ae, be))
end

function check_NaN(X...)
    for x in X
        @show maxabs(x)
        any(isnan, x) && error("NaN encountered")
    end
    error("just terminating")
end