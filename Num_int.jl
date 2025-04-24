using FundamentalsNumericalComputation

μ = 398600.0 #Gravitational parameter of Earth (km^3/s^2)

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

function euler_orbit(y0, tspan) #OPTIONAL euler calculation(if we want extra data points)
    dt = step(tspan)                
    y = copy(y0)
    trajectory = [y0]

    for _ in 2:length(tspan) #skips first point since that is the IC
        dy = rhs(y)
        y = y .+ dt .* dy   #eulers
        push!(trajectory, copy(y))
    end

    return hcat(trajectory...)
end

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


#Initial state: x, y, z, vx, vy, vz
y0 = [20000, -105000, -19000, 0.9, -3.4, -1.5] #TEST CASE FROM TRISTAN

tspan = range(0, stop=60*60*2, length=10000) #Time span: ti, time(seconds), number of iterations


#traj = euler_orbit(y0, μ, tspan)
traj = simpson_orbit(y0, tspan)
display(traj[:, end])

plot3d(traj[1, :], traj[2, :], traj[3, :])