using LinearAlgebra

function cubspline(t,y)
    n = length(t)-1
    h = zeros(1,n)
    for i = 1:n
        h[i] = t[i+1]-t[i]
    end
    Id = I(n)
    z = zeros(n,n)
    E = Id[1:n-1,:]
    J = diagm(0=>ones(n),1=>-ones(n-1))
    first_mat = [Id z z z]
    H = diagm(0=>h)
    second_mat = [Id H H.^2 H.^3]
    third_mat = E*[z J 2*H 3*H.^2]
    fourth_mat = E*[z z J 3*H]
    fifth_mat = [ zeros(1,3*n) [1 -1 zeros(1,n-2)] ]
    sixth_mat = [ zeros(1,3*n) [zeros(1,n-2) 1 -1] ]
    RHS = [y[1:n]; y[2:n+1]; zeros(1,n-1); 0; 0]
    full_mat = [first_mat;second_mat;third_mat;fourth_mat;fifth_mat;sixth_mat]
    coeff = full_mat/RHS

    a = coeff[1:n]
    b = coeff[n+1:2*n]
    c = coeff[2*n+1:3*n]
    d = coeff[3*n+1:end]

    F = [ Polynomial([a[k],b[k],c[k],d[k]]) for k in 1:n ]

    return function (x)
        if x < t[1] || x > t[n+1]    # outside the interval
            return NaN
        elseif x==t[1]
            return y[1]
        else
            k = findlast(x .> t)    # last node to the left of x
            return S[k](x-t[k])
        end
    end

end