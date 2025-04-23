
using LinearAlgebra, Polynomials

function cubspline(t,y)
    n = length(t)-1
    h = zeros(n)
    
    for i = 1:n
        h[i] = t[i+1]-t[i]
    end
    @show size(h)
    Id = I(n)
    z = zeros(n,n)
    E = Id[1:n-1,:]
    J = diagm(0=>ones(n),1=>-ones(n-1))
    first_mat = [Id z z z]
    H = diagm(0 => h)
    @show size(H)
    second_mat = [Id H H.^2 H.^3]
    third_mat = E*[z J 2*H 3*H.^2]
    fourth_mat = E*[z z J 3*H]
    fifth_mat = hcat(zeros(1, 3n), [1.0], [-1.0], zeros(1, n-2))
    sixth_mat = hcat(zeros(1, 3n), zeros(1, n-2), [1.0], [-1.0])    
    @show size(fifth_mat)
    RHS = vcat(y[1:n], y[2:n+1], zeros(2*(n-1)), 0, 0)
    @show size(RHS)
    full_mat = [first_mat;second_mat;third_mat;fourth_mat;fifth_mat;sixth_mat]
    coeff = full_mat\RHS

    a = coeff[1:n]
    b = coeff[n+1:2*n]
    c = coeff[2*n+1:3*n]
    d = coeff[3*n+1:end]

    F = [ Polynomial([a[k],b[k],c[k],d[k]]) for k in 1:n ]

    return function (x)
        if x < t[1] || x > t[n+1]
            return NaN
        end
    
        k = findlast(i -> x > t[i], 1:n)
        if isnothing(k)
            return y[1]
        else
            return F[k](x - t[k])
        end
    end    
end

t = [0.0, 1.0, 2.0, 3.0]
y = [0.0, 1.0, 0.0, 1.0]

spline = cubspline(t, y)
spline2 = spinterp(t,y)



println(spline(0.5))   # Should return an interpolated value between 0 and 1
println(spline(2.5))   # Should return an interpolated value between 0 and 1
println(spline(3.5))   # Outside the domain → returns NaN


# Generate sample points for f(x) = x^2
t = collect(0.0:1.0:9.0)
y = t .^ 2

# Create spline interpolator
spline = cubspline(t, y)

# Evaluate on a dense grid for plotting
x_vals = range(t[1], t[end], length=300)
y_vals = [spline(x) for x in x_vals]

# Plot
plot(x_vals, y_vals, label="Cubic Spline", lw=2, color=:blue)
scatter!(t, y, label="Data Points", color=:red, marker=:circle, ms=4)
xlabel!("x")
ylabel!("f(x)")
title!("Cubic Spline Interpolation of x^2")