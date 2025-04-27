
using FundamentalsNumericalComputation#,plots, DifferentialEquations, LinearAlgebra, Polynomials

"""
------------------------------------------------------------------------------------------------------------------------------------------------------------
The functions within this section are all of the functions made for this project. These include
    RK4 method
    Numerical Integration using Simpson's Rule
    Finite Difference Method
    Cubic Spline Interpolation
------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

"""
    rk4(ivp,n)

The Runge-Kutta 4th order method to solve the given IVP using `n` time steps.
"""
function rk4(ivp,n)
    #extracting the values from the ivp struct
    ti, tfinal = ivp.tspan
    f = ivp.f
    u0 = ivp.u0
    p = ivp.p
    
    h = (tfinal - ti)/n
    t = [ti + i*h for i in 0:n]

    #initalizing the output array
    u = fill(float(u0), n+1)

    #stepping through the time values and applying the RK4 method
    for i in 1:n
        k1 = f(u[i], p, t[i])
        k2 = f(u[i] + (h/2) * k1, p, t[i] + h/2)
        k3 = f(u[i] + (h/2) * k2, p, t[i] + h/2)
        k4 = f(u[i] + h * k3, p, t[i] + h)
        u[i+1] = u[i] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    end
    return t,u
end


"""
    simpson_orbit(y0, tspan)
Numerical integration using Simpson's rule
"""
function simpson_orbit(y0, tspan) #Numerical integration using simpson's rule
    dt = step(tspan)
    y = copy(y0)
    trajectory = [y0]

    for i in 2:length(tspan)
        k1 = rhs(y)
        y_half = y .+ dt/2 .* k1       #midpoint estimate
        k2 = rhs(y_half)
        y_full = y .+ dt .* k2         #full-step estimate
        k3 = rhs(y_full)

        y = y .+ dt/6 .* (k1 .+ 4 .* k2 .+ k3) #Simpson's approximation
        push!(trajectory, copy(y))
    end

    return hcat(trajectory...)
end

"""
    fdEuler(y, t)
Finite Difference Method
"""
function fdEuler(y, t)
    global μ

    # ∆t
    dt = step(t)
    
    # Position Magnitude
    posMag = sqrt((y[1])^2+(y[2])^2+(y[3])^2)
   
    # Acceleration
    xAcc = -μ * y[1]/posMag^3
    yAcc = -μ * y[2]/posMag^3
    zAcc = -μ * y[3]/posMag^3

    # Position Coordinates
    xPos = [y[1], y[1] + dt * y[4] + .5 * dt^2 * xAcc]
    yPos = [y[2], y[2] + dt * y[5] + .5 * dt^2 * yAcc]
    zPos = [y[3], y[3] + dt * y[6] + .5 * dt^2 * zAcc]

    # Velocity Coordinates
    xVel = [y[4], (xPos[end] - xPos[end - 1])/dt]
    yVel = [y[5], (yPos[end] - yPos[end - 1])/dt]
    zVel = [y[6], (zPos[end] - zPos[end - 1])/dt]

    for _ in 3:length(t)
        # Finding New Position Magnitude
        posMag = sqrt((xPos[end])^2+(yPos[end])^2+(zPos[end])^2)

        # Finding New Acceleration
        xAcc = -μ * xPos[end]/posMag^3
        yAcc = -μ * yPos[end]/posMag^3
        zAcc = -μ * zPos[end]/posMag^3

        # Using 2nd-Order Estimates to Find Next Positions
        push!(xPos, 2*xPos[end] - xPos[end - 1] + dt^2 * xAcc)
        push!(yPos, 2*yPos[end] - yPos[end - 1] + dt^2 * yAcc)
        push!(zPos, 2*zPos[end] - zPos[end - 1] + dt^2 * zAcc)

        # Calculating 1st-Order Derivatives
        push!(xVel,(xPos[end] - xPos[end - 1])/dt)
        push!(yVel,(yPos[end] - yPos[end - 1])/dt)
        push!(zVel,(zPos[end] - zPos[end - 1])/dt)
    end
    return hcat(xPos,yPos,zPos,xVel,yVel,zVel)
end


"""
    cubspline(t,y)
"""
function cubspline(t,y)
    n = length(t)-1
    h = zeros(n)
    
    for i = 1:n
        h[i] = t[i+1]-t[i]
    end
    #@show size(h)
    Id = I(n)
    z = zeros(n,n)
    E = Id[1:n-1,:]
    J = diagm(0=>ones(n),1=>-ones(n-1))
    first_mat = [Id z z z]
    H = diagm(0 => h)
    #@show size(H)
    second_mat = [Id H H.^2 H.^3]
    third_mat = E*[z J 2*H 3*H.^2]
    fourth_mat = E*[z z J 3*H]
    fifth_mat = hcat(zeros(1, 3n), [1.0], [-1.0], zeros(1, n-2))
    sixth_mat = hcat(zeros(1, 3n), zeros(1, n-2), [1.0], [-1.0])    
    #@show size(fifth_mat)
    RHS = vcat(y[1:n], y[2:n+1], zeros(2*(n-1)), 0, 0)
    #@show size(RHS)
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

"""
Right-hand side function for Simpson's rule.
"""
function rhs(y)  #setup for the ODE
    x, y_, z, vx, vy, vz = y    #y_ used since y is already being used
    r3 = (x^2 + y_^2 + z^2)^(3/2)
    return [
        vx,
        vy,
        vz,
        -μ * x / r3,
        -μ * y_ / r3,
        -μ * z / r3
    ]
end


"""
    function f(y, μ)

The equations for Kepler's two-body problem are given by the following equations:
    x'' = -μ * x / r^3
    y'' = -μ * y / r^3
    z'' = -μ * z / r^3

Where μ is the gravitational parameter of the Earth (in km^3/s^2) and r is the distance from the center of the Earth to the satellite. 
For the methods, this function transforms the equations into a system of first-order ODEs.
"""
function f(y, μ, t)
    #Where x = y1, y = y2, z = y3 (positions) and vx = y4, vy = y5, vz = y6 (velocities)
    y1, y2, y3, y4, y5, y6 = y 
    r_magitude = sqrt(y1^2 + y2^2 + y3^2)
    y1_dot = y4
    y2_dot = y5
    y3_dot = y6

    #Time derivaties of the velocities (from Kepler's equation)
    y4_dot = -μ * y1 / r_magitude^3
    y5_dot = -μ * y2 / r_magitude^3
    y6_dot = -μ * y3 / r_magitude^3

    #Return array of derivatives
    return [y1_dot, y2_dot, y3_dot, y4_dot, y5_dot, y6_dot]
end  

"""
------------------------------------------------------------------------------------------------------------------------------------------------------------
Case set up and initialization for the methods. This test case is for a 2 hour simulation of a satellite in orbit around the Earth.
------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
#Gravitational parameter of Earth (km^3/s^2)
μ = 398600.0

#Initial state: x, y, z, vx, vy, vz
y0 = [20000, -105000, -19000, 0.9, -3.4, -1.5] #TEST CASE

n1 = 500
n2 = 100 #number of points for cubic spline interpolation plotting
tspan = range(0, 60*60*2, n1) #Time span: ti, time(seconds), number of iterations

#RK4
rungeKuttaIVP = (f = f, u0 = y0, tspan = (first(tspan), last(tspan)), p = μ)
RK4tvals, RK4uvals = rk4(rungeKuttaIVP, n1)

#Simpson's Rule
simpson_result = simpson_orbit(y0, tspan)

#Finite Difference Euler
fdeuler_result = fdEuler(y0, tspan)

#Extract positions
r_rk4 = [sqrt(u[1]^2 + u[2]^2 + u[3]^2) for u in RK4uvals]
v_rk4 = [sqrt(u[4]^2 + u[5]^2 + u[6]^2) for u in RK4uvals]

r_simpson = sqrt.(simpson_result[1, :].^2 .+ simpson_result[2, :].^2 .+ simpson_result[3, :].^2)
v_simpson = sqrt.(simpson_result[4, :].^2 .+ simpson_result[5, :].^2 .+ simpson_result[6, :].^2)

r_fdeuler = sqrt.(fdeuler_result[:, 1].^2 .+ fdeuler_result[:, 2].^2 .+ fdeuler_result[:, 3].^2)
v_fdeuler = sqrt.(fdeuler_result[:, 4].^2 .+ fdeuler_result[:, 5].^2 .+ fdeuler_result[:, 6].^2)

#Interpolation for each method
spl_r_rk4 = cubspline(RK4tvals, r_rk4)
spl_v_rk4 = cubspline(RK4tvals, v_rk4)

spl_r_simpson = cubspline(collect(tspan), r_simpson)
spl_v_simpson = cubspline(collect(tspan), v_simpson)

spl_r_fdeuler = cubspline(collect(tspan), r_fdeuler)
spl_v_fdeuler = cubspline(collect(tspan), v_fdeuler)

#Different n size for plotting
t_small = range(first(tspan), last(tspan), n2)
r_rk4_smooth = [spl_r_rk4(ti) for ti in t_small]
v_rk4_smooth = [spl_v_rk4(ti) for ti in t_small]
r_simpson_smooth = [spl_r_simpson(ti) for ti in t_small]
v_simpson_smooth = [spl_v_simpson(ti) for ti in t_small]
r_fdeuler_smooth = [spl_r_fdeuler(ti) for ti in t_small]
v_fdeuler_smooth = [spl_v_fdeuler(ti) for ti in t_small]


#Plotting the magnitude of the position vector
p1 = plot(t_small, r_rk4_smooth, label="RK4", lw=2, size=(800,600))
plot!(t_small, r_simpson_smooth, label="Simpson", lw=2, linestyle=:dash, markevery=5)
plot!(t_small, r_fdeuler_smooth, label="Finite Difference", lw=1, linestyle=:dash, markevery=20)
xlabel!("Time (s)")
ylabel!("|r| (km)")
title!("Comparison of Methods for Norm of Position Vector")
display(p1)


#Plotting the magnitude of the velocity vector
p2 = plot(t_small, v_rk4_smooth, label="RK4", lw=2, size=(800,600))
plot!(t_small, v_simpson_smooth, label="Simpson", lw=2, linestyle=:dash, markevery=5)
plot!(t_small, v_fdeuler_smooth, label="Finite Difference", linestyle=:dash, markevery=20)
xlabel!("Time (s)")
ylabel!("|v| (km/s)")
title!("Comparison of Methods for Norm of Velocity Vector")
display(p2)


"""
------------------------------------------------------------------------------------------------------------------------------------------------------------
Error and call time Graphing
------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
#initilalizing the rkIVP
tspan = range(0, 60*60*2, n1) #Time span: ti, time(seconds), number of iterations
rungeKuttaIVP = (f = f, u0 = y0, tspan = (first(tspan), last(tspan)), p = μ)

#RK4 with high N value(closest to exact solution we can get)
display("starting exact calculation")
RK4tvals_exact, RK4uvals_exact = rk4(rungeKuttaIVP, 10^3)
r_rk4_exact = [sqrt(u[1]^2 + u[2]^2 + u[3]^2) for u in RK4uvals_exact]
v_rk4_exact = [sqrt(u[4]^2 + u[5]^2 + u[6]^2) for u in RK4uvals_exact]

#Interpolate the exact using cubic spline to find r and v at arbitrary time t
r_interp = cubspline(RK4tvals_exact, r_rk4_exact)
v_interp = cubspline(RK4tvals_exact, v_rk4_exact)

display("Exact calculation done")

#setting up error loop values
n_loop = @. round(Int, 2^(4:0.5:11))
display(n_loop)
err_r = zeros(length(n_loop), 3)
err_v = zeros(length(n_loop), 3)
time_per_n = zeros(length(n_loop), 3)

#loop through the n values generated
for (k, n) in enumerate(n_loop)
    tspan_n = range(0, 60*60*2, n)


    t_rk4_n = @elapsed RK4tvals_n, RK4uvals_n = rk4(rungeKuttaIVP, n-1)
    r_rk4_n = [sqrt(u[1]^2 + u[2]^2 + u[3]^2) for u in RK4uvals_n]
    v_rk4_n = [sqrt(u[4]^2 + u[5]^2 + u[6]^2) for u in RK4uvals_n]
    
    err_r[k, 1] = norm(r_interp.(tspan_n) - r_rk4_n, Inf)
    err_v[k, 1] = norm(v_interp.(tspan_n) - v_rk4_n, Inf)
    time_per_n[k, 1] = t_rk4_n

    t_simpson_n = @elapsed simpson_result_n = simpson_orbit(y0, tspan_n)
    r_simpson_n = sqrt.(simpson_result_n[1, :].^2 .+ simpson_result_n[2, :].^2 .+ simpson_result_n[3, :].^2)
    v_simpson_n = sqrt.(simpson_result_n[4, :].^2 .+ simpson_result_n[5, :].^2 .+ simpson_result_n[6, :].^2)
    
    err_r[k, 2] = norm(r_interp.(tspan_n) - r_simpson_n, Inf)
    err_v[k, 2] = norm(v_interp.(tspan_n) - v_simpson_n, Inf)
    time_per_n[k, 2] = t_simpson_n

    t_fdeuler_n = @elapsed fdeuler_result_n = fdEuler(y0, tspan_n)
    r_fdeuler_n = sqrt.(fdeuler_result_n[:, 1].^2 .+ fdeuler_result_n[:, 2].^2 .+ fdeuler_result_n[:, 3].^2)
    v_fdeuler_n = sqrt.(fdeuler_result_n[:, 4].^2 .+ fdeuler_result_n[:, 5].^2 .+ fdeuler_result_n[:, 6].^2)
    
    err_r[k, 3] = norm(r_interp.(tspan_n) - r_fdeuler_n, Inf)
    err_v[k, 3] = norm(v_interp.(tspan_n) - v_fdeuler_n, Inf)
    time_per_n[k, 3] = t_fdeuler_n

end

#plotting effect of N on function call time
p_time = plot(n_loop, time_per_n,
    label = [L"RK4" L"Num. Integration: Simpsons" L"Central Finite Difference"],
    title = "Effect of N on Function Call Time",xlabel = "N", ylabel = "Time(sec)")
display(p_time)

#error vs N (r vector)
p_error_r = plot(n_loop, err_r, m = :o, label = [L"RK4" L"Num. Integration: Simpsons" L"Central Finite Difference"])
plot!(p_error_r, n_loop, 10 * 10 * n_loop .^ (-4); 
    l = (:dash, :black), 
    label = "4th order",
    xaxis = (:log10, "n"),
    yaxis = (:log10, "max error"),
    title = "Convergence of Numarical Techniques, r",
    size=(800,600))

plot!(p_error_r, n_loop, 10 * 10 * n_loop .^ (-2); 
    l = (:dash, :gray), 
    label = "2nd order")

display(p_error_r)

#error vs N (r vector)
p_error_v = plot(n_loop, err_r, m = :o, label = [L"RK4" L"Num. Integration: Simpsons" L"Central Finite Difference"])
plot!(p_error_v, n_loop, 10 * 10 * n_loop .^ (-4); 
    l = (:dash, :black), 
    label = "4th order",
    xaxis = (:log10, "n"),
    yaxis = (:log10, "max error"),
    title = "Convergence of Numarical Techniques, v",
    size=(800,600))

plot!(p_error_v, n_loop, 10 * 10 * n_loop .^ (-2); 
    l = (:dash, :gray), 
    label = "2nd order")

display(p_error_v)
