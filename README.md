# cvx-mpc-solarsailing

# A Convex Optimization Approach to Solar Sail Station-Keeping Control in Halo Orbits

## Table of Contents
- [About](#about)
- [Organization](#organization)

## About
We present a convex optimization-based station-keeping algorithm designed for long-term station keeping of unstable halo orbits using a solar sail. 
Our controller determines a sail orientation to minimize deviations from a nominal halo orbit. Traditional methods often linearize the solar-sail propulsion model around nominal angles that define the sail orientation, but this can lead to inaccuracies as the model deviates from the linearization point. Instead, we encode the set of possible thrust vectors generated by the nonlinear solar sail model as the boundary of a convex set, which we then relax to arrive at a convex optimization problem. We demonstrate empirically that this relaxation is tight in most cases (i.e. it produces feasible solutions to the original problem) in realistic simulation examples in the Earth-Moon system, validating the effectiveness of this propulsion-free method for long-term stationkeeping. 


## Organization
This repository contains its own Julia 1.10.1 environment specified by the Project.toml and Manifest.toml files. 

The directories in the project are the following: 
- src: algorithm source code
- examples: contains case studies in the Earth-Moon system
- refs: reference trajectories used for the examples