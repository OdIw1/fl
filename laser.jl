propagate_through!(p::Pulse, M::JonesMatrix) = apply_Jones_matrix!(M, p)

function propagate_through!(p::Pulse, fiber::Fiber)
    u_plotX, u_plotY, U_plotX, U_plotY, n_steps, n_steps_rejected, steps =
        rk4ip_vec!(p, fiber)
    print("rk4ip steps: $n_steps\n")    
end

function propagate_through!(p::Pulse, fout::FileOutput)
    if fout.first_iteration
        fout.UX = similar(p.uX)
        fout.UY = similar(p.uY)
        fout.first_iteration = false
    end

    T = calc_T(p.t)   
    spectrum!(p.uX, fout.UX, p.ifft_plan!, T)
    fout.onlyX || spectrum!(p.uY, fout.UY, p.ifft_plan!, T)

    i = fout.iteration
    postfix = length(fout.postfix) > 0 ? "-" * fout.postfix : ""
    fwrite(fout.outdir, "uX" * postfix, i, p.uX)
    fout.onlyX || fwrite(fout.outdir, "uY" * postfix, i, p.uY)
    fwrite(fout.outdir, "vX" * postfix, i, fout.UX)
    fout.onlyX || fwrite(fout.outdir, "vY" * postfix, i, fout.UY)

    fout.iteration += 1
end

function propagate_through!(p::Pulse, F::RectangularSpectralFilter)
    n = length(p.w)
    dn = integer(n * 0.5(1. - F.bandwidth_ratio))
    mid = div(n, 2)

    p.ifft_plan!(p.uX);                     p.ifft_plan!(p.uY)
    for i = (mid-dn+1):(mid+dn)
        p.uX[i] = 0.im;                     p.uY[i] = 0.im
    end
    p.fft_plan!(p.uX);                      p.fft_plan!(p.uY)
end

function propagate_through!(p::Pulse, F::GaussianSpectralFilter)
    p.ifft_plan!(p.uX);                     p.ifft_plan!(p.uY)
    for i = 1:length(p.t)
        T = exp(-(p.w[i] / F.bandwidth_fr)^2)
        p.uX[i] *= T;                       p.uY[i] *= T
    end
    p.fft_plan!(p.uX);                      p.fft_plan!(p.uY)
end

function propagate_through!(p::Pulse, sa::SaturableAbsorber)
    n = length(p.t)
    for i = 1:n
        P = abs2(p.uX[i]) + abs2(p.uY[i])
        A = sa.modulation_depth / (1 + P/sa.saturation_power)
        T = 1 - A
        p.uX[i] *= T
        p.uY[i] *= T
    end
end

function propagate_through!(p::Pulse, s::SESAM)
    n = length(p.t)
    dt = calc_dt(p.t)
    q = q0 = s.modulation_depth
    t = s.recovery_time
    E = s.saturation_power * t
    
    for i = 1:n
        P = abs2(p.uX[i]) + abs2(p.uY[i])
        dq = (q0 - q)/t - q*P/E
        T = 1 - q - s.unsaturable_loss
        q += dq*dt
        
        p.uX[i] *= T
        p.uY[i] *= T
    end
end

function propagate_through!(p::Pulse, c::Coupler)
    n = length(p.t)
    BLAS.scal!(n, c.transmittance, p.uX, 1)
    BLAS.scal!(n, c.transmittance, p.uY, 1)
end

function propagate_through!(p::Pulse, s::PulseSensor)
    dt = calc_dt(p.t)
    n = length(p.t)
    EX = sqr(BLAS.nrm2(n, p.uX, 1)) * dt;           EY = sqr(BLAS.nrm2(n, p.uY, 1)) * dt
    energy = EX + EY
    print("$(s.name)\t E, nJ: $(energy / 1.e-9) = $EX + $EY\n")
end

function run_laser_scheme!(p::Pulse, laser::LaserScheme, n_iter=1)
    for i = 1:n_iter 
        print_with_color(:green, "iteration $i\n")
        for j = 1:length(laser)
            propagate_through!(p, laser[j])
        end
    end
end
