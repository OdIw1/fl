function fname_tsv(base_fname, i)
    base_fname * "-" * lpad(string(i), 4, '0') * ".tsv"
end

function fwrite(fname, data)
    touch(fname)
    f = open(fname, "w+")
    writedlm(f, data)
    close(f)
end

fwrite(outdir, fname, data) = fwrite(joinpath(outdir, fname), data)

fwrite(outdir, fname, i, data) = fwrite(outdir, fname_tsv(fname, i), data)

function mkpath_today(path)
    ispath(path) || mkpath(path)
    isdir(path) || error("$path already exists and is not a directory")
    path_contents = readdir(path)
    dir_base = string(today())
    
    for i in 0:1000
        dir = dir_base * "-" * lpad(string(i), 4, '0')
        if dir ∉ path_contents
            outdir = joinpath(path, dir)
            mkdir(outdir)
            return outdir
        end
    end
    error("$path seems have thousand of data dirs already, try somewhere else")
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
    [:(fwrite(joinpath($outdir, $(string(v, ".tsv"))), $(fun)($v))) for v in vars]
end

function postfix_vars(postfixes, vars)
    vars_postfixed = Symbol[]
    for postfix in postfixes
        vp = [symbol(string(v, postfix)) for v in vars]
        append!(vars_postfixed, vp)
    end
    return vars_postfixed
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
    vars_postfixed = postfix_vars(["X", "Y"], vars)
    ex = out_expr(outdir, fun, vars_postfixed...)
    return Expr(:block, ex...)
end

macro outv(outdir, vars...)
    vars_postfixed = postfix_vars(["X", "Y"], vars)
    ex = out_expr(outdir, identity, vars_postfixed...)
    return Expr(:block, ex...)
end