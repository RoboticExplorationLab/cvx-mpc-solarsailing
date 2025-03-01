#ODE form to use in DifferentialEquations.jl
function ST_and_dynamics_DFJL_wcontrol!(du, u, p, t)
    
    #original CR3BP units
    #du[1:9] = three_body_prob_dynamics_wcontrol(u[1:9])

    system = p[1]

    #scaled units
    du[1:9] = three_body_prob_dynamics_wcontrol_scaled(system, u[1:9])
    
end

function dynamics_wcontrol_integrate(system, x_0,period)
    
    tspan = (0.0, period)
    prob = ODEProblem(ST_and_dynamics_DFJL_wcontrol!, x_0, tspan, [system])
    sol = solve(prob, TsitPap8(), abstol=1e-12, reltol=1e-12)
    
    return sol
    
end

function integrate_halfperiod(system, x_0, u)
    
    #u will be size 20
    
    #x trajectory not delta x 

    xtraj = zeros(6, 21)

    xtraj[:,1] = x_0

    #go through all the controls
    for i=1:20
        
        #integrate the nonlinear dynamics 1 timstep forward
        x_0_non = [xtraj[:,i]; u[:,i]]

        #has the scaling included
        solution_in = dynamics_wcontrol_integrate(system, x_0_non, Δt)

        #this is the next state
        all_state_next = get_state_wcontrol(solution_in)

        xtraj[:,i+1] = all_state_next[1:6, end]
        
    end
    
    #this is size N. 21. Element 21 is the new 1
    return xtraj
        
end

#integrate only one step
function integrate_dt(system, x_0, u)
        
    #x trajectory not delta x 

    x_0_non = [x_0; u]

    #has the scaling included
    solution_in = dynamics_wcontrol_integrate(system, x_0_non, Δt)

    #println("got solution")

    all_state_next = get_state(solution_in)

    #println("size of all state next: ", size(all_state_next))

    #println(all_state_next)
    
    xnext = all_state_next[1:6, end]

    return xnext

end



function integrate_oneperiod(system, x_0, u)
    
    #u will be size 40 since we only want to apply the controls for one rev
    
    #x trajectory not delta x 

    xtraj = zeros(6, N_period)

    xtraj[:,1] = x_0

    for i=1:N_period-1
        
        #integrate the nonlinear dynamics 1 timstep forward
        x_0_non = [xtraj[:,i]; u[:,i]]

        solution_in = dynamics_wcontrol_integrate(system, x_0_non, Δt)
        #this is the next state
        all_state_next = get_state_wcontrol(solution_in)

        xtraj[:,i+1] = all_state_next[1:6, end]
    end
    
    #this is size N. 40 for this case
    return xtraj
        
end

#ODE for the state transition matrix 
function phi_ode(system, z)
    
    #state vector
    x = z[1:6]
    
    #state transition matrix
    Φ = reshape(z[7:end], (6,6))
    
    #Find the dynamics jacobian at the state vector (CR3BP dynamics)
    #F =  ForwardDiff.jacobian(dx -> three_body_prob_dynamics(dx), x)

    #scaled dynamics
    F =  ForwardDiff.jacobian(dx -> three_body_prob_dynamics_scaled(system, dx), x)
    
    #Differential Equation of State transition Matrix
    Φ_dot = F*Φ
    
    return vec(Φ_dot)
    
end

function ST_and_dynamics(system, z)
    x = z[1:6]
    Φ = z[7:end]
    
    x_dot = three_body_prob_dynamics_scaled(system, x)
    Φ_dot = phi_ode(system, z)
    
    z_dot = [x_dot; Φ_dot]
    
    return z_dot
end

#solve using DifferentialEquations.jl
#dynamics and state transition matrix ode
function ST_and_dynamics_DFJL!(du, u, p, t)
    
    system = p[1]

    du[1:6] = three_body_prob_dynamics_scaled(system, u[1:6])
    du[7:end] = phi_ode(system, u)
    
end

#solve using DifferentialEquations.jl
#just dynamics ode 
function justdynamics_DFJL!(du, u, p, t)
    
    system = p[1]

    du[1:6] = three_body_prob_dynamics_scaled(system, u[1:6])
    
end

#gets the state and state transition matrix for the entire solution

function get_state_stm(solution)
    
    N = size(solution.u)[1]

    all_states = zeros(6, N)

    all_stm = zeros(6,6,N)

    for i=1:N
        all_states[:,i] = solution.u[i][1:6]
        all_stm[:,:,i] = reshape(solution.u[i][7:end], (6,6)) 
    end
    
    #all states and all stm are functions of t
    #solution.t is the time
    return all_states, all_stm, solution.t 
end

#gets the state and transition matrix for the entire solution
#for position velocity and control 
function get_state_wcontrol(solution)
    
    N = size(solution.u)[1]

    all_states = zeros(9, N)

    for i=1:N
        all_states[:,i] = solution.u[i][1:9]
    end
    
    #all states and all stm are functions of t
    #solution.t is the time
    return all_states
end

function get_state(solution)
    
    N = size(solution.u)[1]

    all_states = zeros(6, N)

    for i=1:N
        all_states[:,i] = solution.u[i][1:6]
    end
    
    #all states and all stm are functions of t
    #solution.t is the time
    return all_states
end

#this solved until we cross the xz plane. for the differential corrector it is useful 
function diff_solve(system, x_0, period)
    
    #stop when y=0
    condition(u,t,integrator) = u[2]
    
    affect!(integrator) = terminate!(integrator)
    
    cb = ContinuousCallback(condition, affect!)
    
    tspan = (0.0, period)
    prob = ODEProblem(ST_and_dynamics_DFJL!, x_0, tspan, [system])
    sol = solve(prob, TsitPap8(), abstol = 1e-12, reltol=1e-12, callback=cb)
    
    return sol
    
end

#integrate dynamics and state transition matrix
function dynamics_integrate(system, x_0, period)
    
    tspan = (0.0, period)
    prob = ODEProblem(ST_and_dynamics_DFJL!, x_0, tspan, [system])
    sol = solve(prob, TsitPap8(), abstol=1e-12, reltol=1e-12)
    
    return sol
    
end

#integrate just the dynamics
function just_dynamics_integrate(system, x_0, period)
    
    tspan = (0.0, period)
    prob = ODEProblem(justdynamics_DFJL!, x_0, tspan, [system])
    sol = solve(prob, TsitPap8(), abstol=1e-12, reltol=1e-12)
    
    return sol
    
end


#RK4 integrator with a zero order hold
function RK4_integrator(system, x, u)

    f1 = three_body_prob_dynamics_wcontrol_scaled_RK4(system, x, u)
    f2 = three_body_prob_dynamics_wcontrol_scaled_RK4(system, x+0.5*Δt*f1, u)
    f3 = three_body_prob_dynamics_wcontrol_scaled_RK4(system, x+0.5*Δt*f2, u)
    f4 = three_body_prob_dynamics_wcontrol_scaled_RK4(system, x+ Δt*f3, u)

    xnext = x + (Δt/6.0)*(f1+2*f2+2*f3+f4)

    return xnext

end