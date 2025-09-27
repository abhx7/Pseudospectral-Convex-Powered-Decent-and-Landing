# Pseudospectral Convex Optimization for Powered Descent and Landing

This repository implements and replicates the results of Marco Sagliano’s work on **Pseudospectral Convex Optimization for Powered Descent and Landing (PDL)**.  
The project formulates the rocket landing problem as a convex optimization problem using pseudospectral transcription, enabling efficient trajectory optimization under thrust, mass, and glideslope constraints.  

---

## 📖 Description

The code demonstrates how convex optimization can be applied to powered descent and landing problems for reusable launch vehicles (e.g., Falcon 9), providing fuel-optimal and dynamically feasible trajectories.  

It includes:  
- Pseudospectral transcription of the dynamics  
- Convexification of non-convex thrust magnitude constraints  
- Implementation in MATLAB with cvx
- Replication of trajectories shown in Sagliano’s papers  

---

## 📊 Results

The following plots are generated during simulation and optimization:

- **Position (x, y, z) vs Time**
![3D Trajectory](images/trajectory.png) 
- **Velocity (vx, vy, vz) vs Time**
![3D Trajectory](images/trajectory.png)
- **Throttle vs Time**
![3D Trajectory](images/trajectory.png)
- **Mass consumption over time**
![3D Trajectory](images/trajectory.png)  
- **3D Trajectory**
![3D Trajectory](images/trajectory.png)
- **Trajectory animation**  
![3D Trajectory](images/trajectory.png)

---

## 📚 References

- Sagliano, M. *Pseudospectral Convex Optimization for Powered Descent and Landing* (original research work replicated here).  
- Haikal Fouzi’s blog post: [Convex & Optimal Rocket Landing Guidance Control](https://haikalfouzi.com/2021/01/03/Convex-&-Optimal-Rocket-Landing-Guidance-Control.html) — a clear explanation of the approach.  

---
