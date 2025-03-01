#solve.jl

#take out bias, gammas_k, eclipse_k 


#formulate the problem at every solve
#pass in the initial condition, along with jacobians of the current reference trajectory 
function update_prob(system, x_initial_k, all_Ad_k, all_Bd_k, reference_traj_k, sun_satellite_position_k)
    
    X = Convex.Variable(nx,N_h)
        
    U = Convex.Variable(nu, N_h-1)
    
    #initial state constraint
    cons = Constraint[X[:,1] == x_initial_k]


    #dynamics constraint 
    for k=1:(N_h-1)

        push!(cons, zeros(6)== reference_traj_k[:,k+1] + X[:,k+1] - (dynamics_wcontrol_integrate(system, [reference_traj_k[:,k]; zeros(3)], Δt).u[end][1:6] + all_Ad_k[:,:,k]*X[:,k]+all_Bd_k[:,:,k]*U[1:3,k]))

    end

    for k=1:N_h-1

        #with avg gamma. make sure it is consistent with the translation
        P_k_1 = eig_1*sun_satellite_position_k[:,k]*sun_satellite_position_k[:,k]' + eig_2*(Matrix(1.0*I, 3,3) - sun_satellite_position_k[:,k]*sun_satellite_position_k[:,k]')

        P_k_hermitian = (P_k_1 + P_k_1')/2
        
        #constant translation... #works well
        t_k = (avg_gamma/2)*sun_satellite_position_k[:,k]

        U_translated = U[1:3,k] - t_k

        #relaxed ellipsoid constraint
        sail_constraint = quadform(U_translated, P_k_hermitian) <= 1 #+ U[4,k]

        #greater than or equal to constraint is nonconvex
        
        push!(cons, sail_constraint)

    end

    #add a half space constraint
    for k=1:N_h-1

        push!(cons, dot(U[1:3,k], sun_satellite_position_k[:,k]) >= 0)

    end


    return cons, X, U
        
end



#create the optimization problem

function solve_opt(cons, X, U,  N)
    
    #the cost function (L1 norm on the cost and the tracking error)
    obj = 0
    
    for k=1:N

        #L1 norm (working good)
        obj += norm(X[:,k], 1) 
        
    end

    prob = minimize(obj, cons);

    #solve the problem with clarabel
    Convex.solve!(prob, ()->Clarabel.Optimizer(), silent= true)

    Xm = X.value;
    Um = U.value;
    
    return Xm, Um, prob.status
    
end


#check the ellipsoid constraint 
function constraint_check_f(Uk, sun_sat_pose_k)

    constraint_check_1 = zeros(160)

    for k=1:N_h-1
    
        P_k_1 = eig_1*sun_sat_pose_k[:,k]*sun_sat_pose_k[:,k]' + eig_2*(Matrix(1.0*I, 3,3) - sun_sat_pose_k[:,k]*sun_sat_pose_k[:,k]')
    
        P_k_hermitian = (P_k_1 + P_k_1')/2
    
        t_k = (avg_gamma/2)*sun_sat_pose_k[:,k]
    
        U_translated = Uk[1:3,k] - t_k
    
        constraint_check_1[k] = U_translated'*P_k_hermitian*U_translated
    
    
    end

    return constraint_check_1
    
end