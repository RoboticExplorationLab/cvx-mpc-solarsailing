#solve projection.jl

#all the functions for the projection step that we solve with IPOPT


#updated for u and n variables in the ipopt projection problem
#create indices
function create_idx(n_u,n_n,N)
    # This function creates some useful indexing tools for Z 
    # x_i = Z[idx.x[i]]
    # u_i = Z[idx.u[i]]

    #this has just states and controls
    nz = (N-1) * n_u + (N-1) * n_n # length of Z 

    #range (indices) of where the states are in the Z variable
    u = [(i - 1) * (n_u + n_n) .+ (1 : n_u) for i = 1:(N-1)]

    #range (indices) of where the controls are in the Z variable
    n = [(i - 1) * (n_u + n_n) .+ ((n_u + 1):(n_u + n_n)) for i = 1:(N - 1)]


    # constraint indexing for the (N-1) control constraints when stacked up
    c = [(i - 1) * (n_u) .+ (1 : n_u) for i = 1:(N - 1)]
    nc = (N - 1) * n_u
    
    return (n_u=n_u,n_n=n_n,N=N,nz=nz,nc=nc, u= u,n = n,c = c)

end


#L2 norms
function cost(params, Z)

    idx, N = params.idx, params.N

    ellipsoid_u = params.ellipsoid_u

    J = 0

    for i=1:N-1

        ui = Z[idx.u[i]]

        J += (norm(ui - ellipsoid_u[:,i]))^2 

    end

    return J
    
end

function solar_sail_equality_constraint(params, Z::Vector)

    idx, N = params.idx, params.N

    #this should be a 3x(N-1) matrix
    sun_vectors_normal = params.sun_vectors_normal

    avg_gamma_k = params.avg_gamma

    c = zeros(eltype(Z), idx.nc)

    for i=1:(N-1)

        u_i = Z[idx.u[i]]

        n_i = Z[idx.n[i]]
        
        c[idx.c[i]] = avg_gamma_k*((sun_vectors_normal[:,i]'*n_i)^2)*n_i - u_i

    end

    return c

end

function normal_vec_equality_constraint(params, Z::Vector)

    idx, N = params.idx, params.N

    c = zeros(eltype(Z), (N-1))

    for i=1:N-1

        n_i = Z[idx.n[i]]

        c[i] = n_i'*n_i - 1

    end


    return c

end


#this inequality is less than or equal to zero
function dot_product_inequality_constraint(params, Z::Vector)

    idx, N = params.idx, params.N

    #this should be a 3x(N-1) matrix
    sun_vectors_normal = params.sun_vectors_normal

    c = zeros(eltype(Z), (N-1))

    for i=1:N-1

        n_i = Z[idx.n[i]]
    
        c[i] = -(sun_vectors_normal[:,i]'*n_i)

    end

    return c

end

function equality_constraint(params, Z)

    return [
            solar_sail_equality_constraint(params, Z); 
            normal_vec_equality_constraint(params, Z)
            ]

end

function inequality_constraint(params, Z)

    return dot_product_inequality_constraint(params, Z)

end

function set_initial_state(params)

    idx = params.idx

    ellipsoid_u = params.ellipsoid_u

    initial_u = ellipsoid_u
    
    initial_n = zeros(3, idx.N-1)

    for i=1:(idx.N-1)

        initial_n[:,i] = initial_u[:,i]/norm(initial_u[:,i])

    end

    initial_state = zeros(idx.nz)

    #fill in the initial state

    for i=1:size(idx.u)[1]

        initial_state[idx.u[i]] = initial_u[:,i]
    
    end

    for i=1:size(idx.n)[1]

        initial_state[idx.n[i]] = initial_n[:,i]
    
    end

    return initial_state


end

function set_constraint_bounds(params)

    idx = params.idx

    #primal variable bounds (working)
    x_l = -Inf*ones(idx.nz)
    x_u =  Inf*ones(idx.nz)

     #inequality constraint bounds
    c_l = -Inf*ones(idx.N-1)
    c_u =  zeros(idx.N-1)

    #equality constraints have no bounds in this version of the wrapper...

    return [x_l, x_u, c_l, c_u]

end

#get the solution in matrix form

function get_matrix_solution(Z, params)

    idx = params.idx

    U_sol = [Z[idx.u[i]] for i = 1:(idx.N-1)]
    N_sol = [Z[idx.n[i]] for i = 1:(idx.N-1)]

    U_sol_matrix = zeros(3, idx.N-1)
    N_sol_matrix = zeros(3, idx.N-1)

    for i=1:idx.N-1

        U_sol_matrix[:, i] = U_sol[i]
    
    end

    for i=1:idx.N-1

        N_sol_matrix[:, i] = N_sol[i] 
    
    end

    return U_sol_matrix, N_sol_matrix

end