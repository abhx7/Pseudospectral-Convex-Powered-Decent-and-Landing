# Pseudospectral Convex Optimization for Powered Descent and Landing

This repository implements and replicates the results of Marco Sagliano’s work on **Pseudospectral Convex Optimization for Powered Descent and Landing (PDL)**.  
The project formulates the rocket landing problem on the Martian environment as a convex optimization problem using pseudospectral transcription, enabling efficient trajectory optimization under thrust, mass, and glideslope constraints.  

---

## 📖 Description

The code demonstrates how convex optimization can be applied to powered descent and landing problems for reusable launch vehicles, providing fuel-optimal and dynamically feasible trajectories.  

It includes:  
- Pseudospectral transcription of the dynamics  
- Convexification of non-convex thrust magnitude constraints  
- Implementation in MATLAB with cvx
- Replication of trajectories shown in Sagliano’s papers  

---

## 📊 Results

The following plots are generated during simulation and optimization:

- **Position (x, y, z) vs Time**
- **Velocity (vx, vy, vz) vs Time**
- **Throttle vs Time**  
![Dynamics](https://github.com/abhx7/Pseudospectral-Convex-Powered-Decent-and-Landing/blob/main/Results/cvx-opt_mars_descent_results-lobatto.png)

- **3D Trajectory**
![3D Trajectory](https://github.com/abhx7/Pseudospectral-Convex-Powered-Decent-and-Landing/blob/main/Results/cvx-opt_mars_descent_trajectory-lobatto.png)

  
- **Mass consumption over time and height compared to a non-convex solution**
![](https://github.com/abhx7/Pseudospectral-Convex-Powered-Decent-and-Landing/blob/main/Results/cvx-opt_mars_descent_trajectory-lobatto.png)  
![](https://github.com/abhx7/Pseudospectral-Convex-Powered-Decent-and-Landing/blob/main/Results/opt_mars_descent_mass_compareL_z.png)


- **Trajectory animation** 
![3D Trajectory](https://github.com/abhx7/Pseudospectral-Convex-Powered-Decent-and-Landing/blob/main/Results/cvx_traj_animation-lobatto.gif)

---

## 📚 References

- Sagliano, M. *Pseudospectral Convex Optimization for Powered Descent and Landing* (original research work replicated here).  
- Haikal Fouzi’s blog post: [Convex & Optimal Rocket Landing Guidance Control](https://haikalfouzi.com/2021/01/03/Convex-&-Optimal-Rocket-Landing-Guidance-Control.html) — a clear explanation of the approach.  

---
