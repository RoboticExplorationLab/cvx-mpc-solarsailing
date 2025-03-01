#transformations 

#get the jacobian that maps from rotating frame to inertial frame
function get_transformation(tz)

    #time is already wrt et
    
    #state of the moon (position and velocity) relative to Earth (in km and km/s)
    moon_state = spkezr("moon",tz,"J2000","NONE","earth")[1]

    #position of the moon relative to Earth
    r_moon = moon_state[1:3]

    #velocity of the moon relative to Earth
    v_moon = moon_state[4:6]

    #x component of the rotating frame expressed in inertial coordinates
    x̃ = r_moon/norm(r_moon)
    #z component of the rotating frame expressed in inertial coordinates
    z̃ = cross(r_moon, v_moon)/norm(cross(r_moon, v_moon))
    #x component of the rotating frame expressed in inertial coordinates
    ỹ = cross(z̃, x̃)

    #instantanous rotation matrix from rotating frame to inertial frame (centered at Earth)
    C = [x̃ ỹ z̃]
    
    return C  

end

#here x is just a position
function cr3bp_to_eci(x, t)

    #et does not change
    tz = et + t

    #get transformation takes in time wrt et
    C = get_transformation(tz)

    #apply the rotation to get the state in the inertial frame
    x_eci = C*x

    return x_eci

end

#here x is just a position
function eci_to_cr3bp(system, x, t)

    #et does not change
    tz = et + t

    #get transformation takes in time wrt et. Rotating frame to ECI
    C = get_transformation(tz)

    #apply the rotation to get the state in the rotating frame. It's transposed bc C is from cr3bp to eci 
    x_cr3bp_earth_centered = C'*x

    x_cr3bp_barycentric = x_cr3bp_earth_centered - [system.μ*system.position_scale, 0, 0] 

    return x_cr3bp_barycentric, x_cr3bp_earth_centered

end