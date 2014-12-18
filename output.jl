function fname_tsv(base_fname, i)
    lpad(string(i), 4, '0') * "-" * base_fname * ".tsv"
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
    
    for i in 1:1000
        dir = dir_base * "-" * lpad(string(i), 4, '0')
        if dir ∉ path_contents
            outdir = joinpath(path, dir)
            mkdir(outdir)
            return outdir
        end
    end
    error("$path seems to have thousand of data dirs already, try somewhere else")
end

function fread_complex(fname)
    f = open(fname, "r")
    str_data = readdlm(f, '~', String) # '~' hopefully shouldn't be there
    map!(strip, str_data)
    str_data = filter(s -> !isempty(s), str_data)

    n = length(str_data)
    r = zeros(Complex{Float64}, n)
    for i = 1:n
        cstr = str_data[i]
        cstr_splitted = split(cstr, " ")
        re_str, im_str = first(cstr_splitted), last(cstr_splitted)
        im_str = replace(im_str, "im", "")
        sign = ("+" in cstr_splitted) ? 1.0 : -1.0
        r[i] = Complex(float64(re_str), sign * float64(im_str))
    end
    r
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



function preprocess_plot(dir::String, fname_postfix::String;
                         t_points=0, t_low=0, t_up=0, i_low=0, i_up=0)
    isdir(dir) || error("$dir is not a valid directory")
    dir_contents = readdir(dir)
    
    R = Regex("^.*" * fname_postfix * ".tsv\$")
    files = sort(filter(f -> ismatch(R, f), dir_contents))

    n_files = length(files)
    i_low = i_low == 0 ? 1: clamp(i_low, 1, n_files)
    i_up = i_up == 0 ? n_files: clamp(i_up, 1, n_files) 
    files = files[i_low:i_up]

    processed = {}
    for fname in files
        print("processing $fname\r") 
        d = fread_complex(joinpath(dir, fname))
        r = resample(d, t_points, t_low, t_up)
        print(r)
        push!(processed, r)
    end

    outname = "!" * fname_postfix * "_resampled.tsv"
    f = open(joinpath(dir, outname), "w+")
    writedlm(f, processed , ',')
    close(f)
    
    return processed
end

function resample{T}(a::Vector{T}, t_points=0::Integer, t_low=0::Integer, t_up=0::Integer)
    n = length(a)
    t_low    = t_low == 0 ? 1: clamp(t_low, 1, n)
    t_up     = t_up == 0 ? n: clamp(t_up, 1, n)
    if t_points == 0 
        return a
    else 
        t_points = clamp(t_points, 1, n)
    end

    averaging_interval_width = int(ceil((t_up - t_low)/(2*(t_points-1))))
    t_indices = map(int, round(linspace(t_low, t_up, t_points)))
    r = similar(t_indices, T)
    
    for i = 1:t_points
        l = max(1, t_indices[i] - averaging_interval_width)
        u = min(n, t_indices[i] + averaging_interval_width)
        r[i] = mean(a[l:u])
    end
    return r
end 
