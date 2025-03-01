#Resources to obtain the position scale parameters 
# https://ross.aoe.vt.edu/books/Ross_3BodyProblem_Book_2022.pdf
# https://ssd.jpl.nasa.gov/tools/periodic_orbits.html


#define a struct for the three body system
struct ThreeBodySystem

    #mass paramter
    μ::Float64

    #constants dependent on μ
    #μ1 = 1-μ
    μ1::Float64
    #μ2 = μ
    μ2::Float64

    mass::Float64

    #position of mass 1
    pose_m1::Vector{Float64}

    #position of mass 2
    pose_m2::Vector{Float64}

    #position scale
    position_scale::Float64

    #timescale 
    time_scale::Float64

    #velocity scale
    velocity_scale::Float64

    #acceleration scale 
    acceleration_scale::Float64

    #position of L1 Lagrange Point
    XL1::Vector{Float64}

    #position of L2 Lagrange Point 
    XL2::Vector{Float64}

end


#create a custom constructor Earth-Moon to define all the parameters and save them in a struct
function ThreeBodySystem_EarthMoon()

    μ = 1.215e-2

    μ1 = 1-μ

    μ2 = μ

    #mass of the spacecraft (kg)
    mass = 4 

    pose_m1 = [-μ, 0.0, 0.0]

    pose_m2 = [1-μ, 0.0, 0.0]

    #these values are from the following textbook linked above
    #positions will be in km 
    position_scale = 3.850e5

    #time will be in days 
    time_scale = (2.361e6/(86400))/(2*pi)

    velocity_scale = position_scale/time_scale

    acceleration_scale = position_scale/(time_scale^2)

    #these are computed offline
    XL1 = [0.8369180073169304, 0, 0, 0, 0, 0]

    XL2 = [1.1556799130947355, 0, 0, 0, 0, 0]

    return ThreeBodySystem(μ, μ1, μ2, mass, pose_m1, pose_m2, position_scale, time_scale, velocity_scale, acceleration_scale, XL1, XL2)

end

