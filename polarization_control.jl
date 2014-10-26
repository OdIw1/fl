const half_wave = [-1.im 0; 0 1.im]
const quarter_wave = 1. / sqrt(2) * [(1. - 1.im)  0; 0 (1. + 1.im)]

function R(a)
	[cos(a)  -sin(a); sin(a) -cos(a)]
end

function W2(a)
	R(a) * half_wave * R(-a)
end

function W4(a)
	R(a) * quarter_wave * R(-a)
end

function apply_Jones_matrix!(M, uX, uY)
    for i = 1:lenght(uX)
        uX_ = uX[i]                 
        uY_ = uY[i]
        uX[i] = a[1,1] * uX_ + a[1,2] * uY_
        uY[i] = a[2,1] * uX_ + a[2,2] * uY_
    end
end