function ds = mars_descent_dynamics(t, s, U, g_vec, Tmax, n_thrusters, phi, alpha)
% mars_descent_dynamics
% Computes the time derivative of the state during Mars powered descent.
%
% Inputs:
%   t       - Time (unused but required for ODE compatibility)
%   s       - State vector [r; v; m] ∈ ℝ⁷, where:
%               r ∈ ℝ³ = position [rx; ry; rz] in meters
%               v ∈ ℝ³ = velocity [vx; vy; vz] in m/s
%               m ∈ ℝ  = mass in kg
%   U       - Normalized thrust control vector [Tx; Ty; Tz], each ∈ [Tl, Tu]
%   g_vec   - Constant gravity vector, e.g., [0; 0; -3.7114] m/s²
%   Tmax    - Maximum thrust per thruster [N]
%   phi     - Cant angle of each thruster [rad]
%   alpha   - Mass flow coefficient: α = 1 / (Isp * gₑ * cos(phi))
%
% Output:
%   ds      - State derivative [dr; dv; dm] ∈ ℝ⁷

    % Extract state variables
    r = s(1:3);             % Position vector [m]
    v = s(4:6);             % Velocity vector [m/s]
    m = s(7);               % Mass [kg]

    % Compute thrust components from normalized control U
    Tcx = Tmax * n_thrusters * cos(phi) * U(1);  % Thrust along x
    Tcy = Tmax * n_thrusters * cos(phi) * U(2);  % Thrust along y
    Tcz = Tmax * n_thrusters * cos(phi) * U(3);  % Thrust along z

    Tc = [Tcx; Tcy; Tcz];    % Total thrust vector [N]

    % Translational dynamics
    drdt = v;                          % ṙ = v
    dvdt = Tc / m + g_vec;            % v̇ = T/m + g

    % Mass depletion rate
    dmdt = -alpha * norm(Tc);         % ṁ = -α * ||T||

    % Return state derivative
    ds = [drdt; dvdt; dmdt];
end
