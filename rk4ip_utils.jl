type RK4IPTemp{T<:Real}
    n::Integer 
    k1X::Vector{Complex{T}}
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
                                              ue_cplx_::Vector{T}, atol=1.e-6, rtol=1.e-6)
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