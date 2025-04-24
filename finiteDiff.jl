using FundamentalsNumericalComputation

function fdEuler(y, t)
    # Gravitational parameter of Earth (km^3/s^2)
    μ = 398600.0

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


# Initial state: x, y, z, vx, vy, vz
y0 = [20000, -105000, -19000, 0.9, -3.4, -1.5]

# Time span: ti, time(seconds), number of iterations
tspan = range(0, stop=60*60*2, length=10000)

traj = fdEuler(y0, tspan)
display(traj[end,:])

plot3d(traj[:,1], traj[:,2], traj[:,3])