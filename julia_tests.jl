function test_argument_passing()
    x = [1, 2]
    f1(x)
    x
end

function f1(x)
    x = 100
end

function run_closures()
    test_closure_bindings()();
    test_closure_bindings2()();
    test_closure_bindings3()();       
end

function f_test(x, y)
    println("x: $(x), y: $(y)")
end

function create_closure(x, y)
    f = () -> f_test(x, y)
end

function test_closure_bindings()
    x = y = [1, 2, 3]    
    f = () -> f_test(x, y)
    x = -10
    y[1] = -10
    f
end

function test_closure_bindings2()
    x = y = [1, 2, 3]
    f = create_closure(x, y)
    x =  -10
    y[1] = -10
    f
end

function test_closure_bindings3()
    x = y = [1, 2, 3]
    f = cell(3)
    let x = x, y = y  
        f[1] = () -> f_test(x, y)
    end
    x[1] = -10
    f[1]
end

function test_fft(k = 6000, n = 2^12)
    u0 = randn(n) + 1.im * randn(n)
    u = copy(u0)
    u_inplace = copy(u0)

    u_plan = copy(u)
    
    FFTW.set_num_threads(4)
    fft_plan = plan_fft(u_plan, (1,), FFTW.ESTIMATE)
    ifft_plan = plan_ifft(u_plan, (1,), FFTW.ESTIMATE)

    # in-place transforms should be faster
    fft_plan! = plan_fft!(u_plan, (1,), FFTW.ESTIMATE)
    ifft_plan! = plan_ifft!(u_plan, (1,), FFTW.ESTIMATE)   

    @profile begin
        tic()
        @time for i in 1:k
            u = fft_plan(ifft_plan(u))
        end
        t = toc()
        @show maximum(abs2(u - u0))

        tic()
        @time for i in 1:k
            ifft_plan!(u_inplace)
            fft_plan!(u_inplace)        
        end
        t! = toc()

        @show maximum(abs2(u_inplace - u0))
        @show inplace_speedup = t / t!
    end
end

function test_slab()
    args = (1, 2)
    f_slab(args...)
end

function f_slab(x,y)
    x + y
end


