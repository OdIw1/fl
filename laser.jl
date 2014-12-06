propagate_through!(p::Pulse, M::JonesMatrix) = apply_Jones_matrix!(M, p)

function propagate_through!(p::Pulse, fiber::Fiber)
    u_plotX, u_plotY, U_plotX, U_plotY, n_steps, n_steps_rejected, steps =
        rk4ip_vec!(p, fiber)    
end

function propagate_through!(p::Pulse, fout::FileOutput)
    UX = similar(p.uX)
    UY = similar(p.uY)
    T = calc_T(p.t)   
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
