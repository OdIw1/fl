propagate_through!(p::Pulse, M::JonesMatrix) = apply_Jones_matrix!(M, p)

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
    postfix = "-" * fout.postfix
    fwrite(outdir, "uX" * postfix, i, p.uX)
    fwrite(outdir, "uY" * postfix, i, p.uY)
    fwrite(outdir, "UX" * postfix, i, UX)
    fwrite(outdir, "UY" * postfix, i, UY)

    fout.iteration += 1
end

function run_laser_scheme!(p::Pulse, laser::LaserScheme, n_iter=2)
    for i = 1:n_iter, j = 1:length(laser)
        propagate_through!(p, laser[j])
    end
end

function run_laser(n, T_window, scheme, L, T0, P0, C0, theta, shape=0)
    scheme = LaserElement[]

    outdir = mkpath_today()
    fout1 = FileOutput(outdir)

    PC1 = QuarterWavePlate(a1)
    _P2 = HalfPlate(a2)
    _P3 = QuarterWavePlate(a3)
    _P4 = Polarizer()
    PC2 = _P4 * _P3 * _P2

    # fiber parameters are from 13[Yarutkina, Shtyrina]{Opt.Expr} Numerical Modeling of ...
    #fiber_active = Fiber(2., 0., beta_fsmm_to_sm([76.9, 168.], 0., 9.32e-3,    )) 
    fiber_passive = FIber(100., 10^0.02 * 1.e-3, beta_fsmm_to_sm([4.5, 109]), 2.1e-3)


end

