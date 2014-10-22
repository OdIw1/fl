function fwrite(fname, data)
    touch(fname)
    f = open(fname, "w+")
    writedlm(f, data)
    close(f)
end

function clampl{T1<:Complex, T2<:Real}(x::T1, lo::T2)
    ax = abs(x)
    ax < lo ? lo : ax
end

function clampl{T<:Complex}(x::T, lo::T)
    ax = abs2(x)
    alo = abs2(lo)
    ax < alo ? alo : ax
end

function clampl{T<:Real}(x::T, lo::T)
    x < lo ? lo : x
end

type ClamplFun <: Functor{2} end
NumericExtensions.evaluate(::ClamplFun, x, y) = clampl(x, y)

function clamp_plot(u, threshold=1.e-4)
    uabs = abs2(u)
    m = threshold * maximum(uabs)
    map1!(ClamplFun(), uabs, m)
    return uabs  
end

function clamp_log_plot(u, threshold=1.e-4)
    uc = clamp_plot(u, threshold)
    map1!(Log10Fun(), uc)
    return uc   
end

function out_expr(outdir, fun, vars...)
    ex = [:(fwrite(joinpath($outdir, $(string(v, ".tsv"))), $(fun)($v))) for v in vars]
    return ex
end

macro outf(outdir, fun, vars...)
    ex = out_expr(outdir, fun, vars...)
    return Expr(:block, ex...)
end

macro out(outdir, vars...)
    ex = out_expr(outdir, identity, vars...)
    return Expr(:block, ex...)
end 

macro outfv(outdir, fun, vars...)
    vars_postfixed = Symbol[]
    for postfix in ["X", "Y"]
        vp = [symbol(string(v, postfix)) for v in vars]
        append!(vars_postfixed, vp)
    end
    ex = out_expr(outdir, fun, vars_postfixed...)
    return Expr(:block, ex...)
end

macro outv(outdir, vars...)
    vars_postfixed = Symbol[]
    for postfix in ["X", "Y"]
        vp = [symbol(string(v, postfix)) for v in vars]
        append!(vars_postfixed, vp)
    end
    ex = out_expr(outdir, identity, vars_postfixed...)
    return Expr(:block, ex...)
end