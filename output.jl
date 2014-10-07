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

function clampl(x, lo)
    x < lo ? lo : x
end

type ClamplFun <: Functor{2} end
NumericExtensions.evaluate(::ClamplFun, x, y) = clampl(x, y)

function clamp_plot(u, threshold=1.e-4)
    uabs = abs(u)
    m = threshold * maximum(uabs)
    map1!(ClamplFun(), uabs, m)
    return uabs  
end

function clamp_log_plot(u, threshold=1.e-4)
    uc = clamp_plot(u, threshold)
    map1!(Log10Fun(), uc)
    return uc   
end