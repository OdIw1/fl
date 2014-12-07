function wl2fr(wl) 2pi * 3.e8 / wl end  
# inverse calculation is the same
function fr2wl(fr) 2pi * 3.e8 / fr end

function bandwidth_wl2fr(wl, dwl) abs(wl2fr(wl + dwl/2) - wl2fr(wl - dwl/2)) end

function bandwidth_wl(n, T, wl0)
    dfr = pi * n / (2T)
    [fr2wl(wl2fr(wl0) + dfr), fr2wl(wl2fr(wl0) - dfr)]
end

function betha_convert(betha, T_factor, L_factor)
    [betha[k] * T_factor^(k+1) / L_factor for k = 1:length(betha)]
end    

ps_km2s_m(betha) = betha_convert(betha, 1.e-12, 1.0e3)
ps_m2s_m(betha)  = betha_convert(betha, 1.e-12, 1.0)
fs_mm2s_m(betha) = betha_convert(betha, 1.e-15, 1.0e-3)

function pulse_params(T0, P0, betha::Array, gamma)
    ld = zeros(Float64, length(betha))
    soliton_order = zeros(ld)

    ln = 1. / (gamma * P0)
    for k = 1:length(betha)
        ld[k] = T0^(k+1) / abs(betha[k])
        soliton_order[k] = sqrt(ld[k] / ln)
    end
    ld, ln, soliton_order
end

function t_grid(n::Integer, T::Real)
# In MATLAB version 'T' was a scaling factor: time window was (-T0 T, T0 T),
# where T0 is pulse FWHM;  in Julia time window is (-T, T)
# also 'n' is no longer binary log of the number of grid points, it is the number itself
    dt = 2 * T / n
    dt * [-n/2 : (n/2 - 1)]
end

function w_grid(n::Integer, T::Real)
    dw = pi / T
    dw * [0 : (n/2 -1), -n/2 : -1]
end

function calc_dt(t)
    n = length(t)
    dt = (t[end] - t[1]) / (n - 1)
end

function calc_T(t)
    dt = calc_dt(t)
    T = (dt + t[end] - t[1]) / 2
end

function secant_pulse(T0, P0, C0, t::Vector, t_offset=0)
    t_scaled = (t - t_offset) / T0
    sqrt(P0) * sech(t_scaled) .* exp(-0.5im * C0 * (t_scaled).^2)
end

function gaussian_pulse(m_order, T0, P0, C0, t::Vector, t_offset=0)
    t_scaled = (t - t_offset) / T0
    sqrt(P0) * exp(-0.5(1 + 1im * C0) * t_scaled .^ (2m_order))
end

function pulse(m_order, T0, P0, C0, t::Vector, t_offset=0)
    if m_order == 0
        secant_pulse(T0, P0, C0, t, t_offset)
    else
        gaussian_pulse(m_order, T0, P0, C0, t, t_offset)
    end
end

function pulse_vec(m_order, T0, P0, C0, theta, t::Vector, t_offset=0)
    u = pulse(m_order, T0, P0, C0, t, t_offset)
    (u * cos(theta), u * sin(theta))
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

macro closure1(name, fun, closure_arg, args...)
    if closure_arg ∉ args
        args_string = join([string(a) for a in args], ",")
        msg = string("closure arg ", closure_arg, " not in ", args_string)
        :(ArgumentError($msg))
    else
        quote
            # global $name
            $(esc(name)) = let $([:($(esc(a)) = $a) for a in args]...)
                ($closure_arg)->$fun($(args...))
            end
        end
    end
end