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
    postfix = length(fout.postfix) > 0 ? "-" * fout.postfix : ""
    fwrite(fout.outdir, "uX" * postfix, i, p.uX)
    fwrite(fout.outdir, "uY" * postfix, i, p.uY)
    fwrite(fout.outdir, "vX" * postfix, i, UX)
    fwrite(fout.outdir, "vY" * postfix, i, UY)

    fout.iteration += 1
end

function run_laser_scheme!(p::Pulse, laser::LaserScheme, n_iter=1)
    for i = 1:n_iter, j = 1:length(laser)
        propagate_through!(p, laser[j])
    end
end

function Yarutkina13(a1=pi/4, a2=0, a3=0)
    laser = LaserElement[]

    PC1 = QuarterWavePlate(a1)
    push!(laser, PC1)

    # fiber parameters are from 13[Yarutkina, Shtyrina]{Opt.Expr} Numerical Modeling of ...
    wl = 1550e-9
    gain_bw_wl = 50e-9
    gain_bw = bandwidth_wl2fr(wl, gain_bw_wl)
    
    t_round = 102. * 1.47 / 3.e8
    sat_e = 20.e-3 * t_round

    fiber_active = Fiber(2., 0., betha_fsmm_to_sm([76.9, 168.]), 9.32e-3, 10^(5.4/10), gain_bw, sat_e)
    push!(laser, fiber_active)

    fiber_passive = Fiber(100., 10^(0.2/10) * 1.e-3, betha_fsmm_to_sm([4.5, 109]), 2.1e-3)
    push!(laser, fiber_passive)

    outdir = mkpath_today("/mnt/hgfs/VM_shared/out")
    fout1 = FileOutput(outdir)
    push!(laser, fout1)

    _P2 = HalfWavePlate(a2)
    _P3 = QuarterWavePlate(a3)
    _P4 = Polarizer()
    PC2 = _P4 * _P3 * _P2
    push!(laser, PC2)

    # seed pulse params
    n = 2^14
    T0 = 1.e-12
    T = 200T0
    t = t_grid(n, T)
    w = t_grid(n, T)
    p = Pulse(T0, 1.e-2, 0., 0., t, w)

    # save the initial pulse
    fout1.iteration = 0
    propagate_through!(p, fout1)

    # run
    run_laser_scheme!(p, laser, 10)
end
