
using Plots, FundamentalsNumericalComputation, DifferentialEquations, LinearAlgebra, Polynomials

"""
------------------------------------------------------------------------------------------------------------------------------------------------------------
All of the functions in this section are from the textbook and are used for comparison purposes.
------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

"""
    spinterp(t,y)

Construct a cubic not-a-knot spline interpolating function for data
values in `y` given at nodes in `t`.
"""
function spinterp(t,y)
    n = length(t)-1
    h = [ t[k+1]-t[k] for k in 1:n ]

    # Preliminary definitions.
    Z = zeros(n,n);
    In = I(n);  E = In[1:n-1,:];
    J = diagm(0=>ones(n),1=>-ones(n-1))
    H = diagm(0=>h)

    # Left endpoint interpolation:
    AL = [ In Z Z Z ]
    vL = y[1:n]

    # Right endpoint interpolation:
    AR = [ In H H^2 H^3 ];
    vR = y[2:n+1]

    # Continuity of first derivative:
    A1 = E*[ Z J 2*H 3*H^2 ]
    v1 = zeros(n-1)

    # Continuity of second derivative:
    A2 = E*[ Z Z J 3*H ]
    v2 = zeros(n-1)

    # Not-a-knot conditions:
    nakL = [ zeros(1,3*n) [1 -1 zeros(1,n-2)] ]
    nakR = [ zeros(1,3*n) [zeros(1,n-2) 1 -1] ]

    # Assemble and solve the full system.
    A = [ AL; AR; A1; A2; nakL; nakR ]
    v = [ vL; vR; v1; v2; 0; 0 ]
    z = A\v

    # Break the coefficients into separate vectors.
    rows = 1:n
    a = z[rows]
    b = z[n.+rows];  c = z[2*n.+rows];  d = z[3*n.+rows]
    S = [ Polynomial([a[k],b[k],c[k],d[k]]) for k in 1:n ]

    # This function evaluates the spline when called with a value
    # for x.
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

Apply the common Runge-Kutta 4th order method to solve the given 
IVP using `n` time steps. Returns a vector of times and a vector of
solution values.
"""
function rk4(ivp,n)
    # Time discretization.
    a,b = ivp.tspan
    h = (b-a)/n
    t = [ a + i*h for i in 0:n ]

    # Initialize output.
    u = fill(float(ivp.u0),n+1)

    # Time stepping.
    for i in 1:n
        k₁ = h*ivp.f( u[i],      ivp.p, t[i]     )
        k₂ = h*ivp.f( u[i]+k₁/2, ivp.p, t[i]+h/2 )
        k₃ = h*ivp.f( u[i]+k₂/2, ivp.p, t[i]+h/2 )
        k₄ = h*ivp.f( u[i]+k₃,   ivp.p, t[i]+h   )
        u[i+1] = u[i] + (k₁ + 2(k₂+k₃) + k₄)/6
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

    for _ in 2:length(tspan)
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
spl_r_rk4 = spinterp(RK4tvals, r_rk4)
spl_v_rk4 = spinterp(RK4tvals, v_rk4)

spl_r_simpson = spinterp(collect(tspan), r_simpson)
spl_v_simpson = spinterp(collect(tspan), v_simpson)

spl_r_fdeuler = spinterp(collect(tspan), r_fdeuler)
spl_v_fdeuler = spinterp(collect(tspan), v_fdeuler)

#Different n size for plotting
t_small = range(first(tspan), last(tspan), n2)
r_rk4_smooth = [spl_r_rk4(ti) for ti in t_small]
v_rk4_smooth = [spl_v_rk4(ti) for ti in t_small]
r_simpson_smooth = [spl_r_simpson(ti) for ti in t_small]
v_simpson_smooth = [spl_v_simpson(ti) for ti in t_small]
r_fdeuler_smooth = [spl_r_fdeuler(ti) for ti in t_small]
v_fdeuler_smooth = [spl_v_fdeuler(ti) for ti in t_small]

#Plotting the magnitude of the position vector
p1 = plot(t_small, r_rk4_smooth, label="RK4", lw=2, color=:blue)
plot!(t_small, r_simpson_smooth, label="Simpson", lw=2, color=:red, linestyle=:dash)
plot!(t_small, r_fdeuler_smooth, label="FD Euler", lw=2, color=:green, linestyle=:dash)
#scatter!(RK4tvals, r_rk4, label="RK4 Points", marker=:circle, color=:blue)
#scatter!(collect(tspan), r_simpson, label="Simpson Points", marker=:square, color=:red)
#scatter!(collect(tspan), r_fdeuler, label="FD Euler Points", marker=:diamond, color=:green)
xlabel!("Time (s)")
ylabel!("|r| (km)")
title!("Comparison of Methods for Norm of Position Vector")
display(p1)


#Plotting the magnitude of the velocity vector
p2 = plot(t_small, v_rk4_smooth, label="RK4", lw=2, color=:blue)
plot!(t_small, v_simpson_smooth, label="Simpson", lw=2, color=:red, linestyle=:dash)
plot!(t_small, v_fdeuler_smooth, label="FD Euler", lw=2, color=:green, linestyle=:dash)
#scatter!(RK4tvals, v_rk4, label="RK4 Points", marker=:circle, color=:blue)
#scatter!(collect(tspan), v_simpson, label="Simpson Points", marker=:square, color=:red)
#scatter!(collect(tspan), v_fdeuler, label="FD Euler Points", marker=:diamond, color=:green)
xlabel!("Time (s)")
ylabel!("|v| (km/s)")
title!("Comparison of Methods for Norm of Velocity Vector")
display(p2)








