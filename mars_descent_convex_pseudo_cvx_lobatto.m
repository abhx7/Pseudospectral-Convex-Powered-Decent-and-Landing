%% mars_pdl_cvx.m
% CVX-based solve of the convexified Mars powered descent problem using pseudospectral transcription.

clear; clc; close all;

%% ------------------ Parameters  ------------------
g        = [0; 0; -3.7114];            % Mars gravity [m/s^2]
Isp      = 225;                        % Specific impulse [s]
ge       = 9.807;                      % Earth gravity [m/s^2]
phi      = deg2rad(27);                % Thruster cant [rad]
Tmax_ax  = 3100;                       % Max thrust per axis per thruster [N]
Tl       = 0.3;  Tu = 0.8;             % Throttle bounds (fraction of Tmax)
thetaAlt = deg2rad(4);                 % Glideslope angle [rad]
nThr     = 6;                          % Number of thrusters

alpha    = 1/(Isp*ge*cos(phi));        % Fuel flow coeff 
TmaxTot  = nThr*Tmax_ax*cos(phi);               % cap on ||T|| (vector magnitude)

ns = 7;   % number states per node (rx,ry,rz,vx,vy,vz,z)
nc = 4;   % controls per node (ux,uy,uz,sigma)

m_dry = 1505;

% Time grid
t0 = 0;
tf = 81;                               % Final time [s]
N  = 30;                               % Number of collocation nodes (decision nodes)

% Initial / final conditions
r0 = [2000; 0; 1500];
v0 = [100;  0; -75];
m0 = 1905;

rf = [0; 0; 0];
vf = [0; 0; 0];

% % Simple initial guess for continuity/plotting (not used by CVX solver, but useful)
% U0 = [-0.32*TmaxTot; 0; 0.48*TmaxTot];
%
% %% ------------------ optional: compute an ODE rollout to visualise initial guess (copied) ------------------
% dynODE = @(t_, s_) [ ...
%     s_(4:6); ...
%     U0./s_(7) + g; ...
%     -alpha*norm(U0) ];
% 
% [t_ode, s_ode] = ode45(dynODE, [t0 tf], [r0; v0; m0]);
% % Interpolate to a uniform pseudospectral-like grid later if needed.

%% ------------------ Pseudospectral nodes and differentiation ------------------
% Use lobatto_nodes helper that returns nodes on [-1,1] and differentiation matrix.
tau = double(lobatto_nodes(N))';   % length N
%tau = double(CGL_Roots(N))';

% Map tau ∈ [-1,1] -> t ∈ [t0,tf]
kt = (tf - t0) / 2;
t_nodes = kt * (tau + 1);   % collocation nodes in [0,tf]

% differentiate matrix
%d = DifferentiationMatrix(N,tau);        % (N)x(N) barycentric diff matrix
D = diff_mat_barycentric(tau);

%% ------------------ Build linear dynamics constraint A0*Y = b0  ------------------
% System matrices for the local ODE: xdot = A x + B u + b_g

Ac = [zeros(3,3), eye(3), zeros(3,1);
      zeros(3,3), zeros(3,3), zeros(3,1);
      zeros(1,3), zeros(1,3), 0];    % ns x ns
Bc = [zeros(3,3), zeros(3,1);
      eye(3), zeros(3,1);
      zeros(1,3), -alpha];          % ns x nc
% Compose block rows to build the big A0 matrix (ns*N rows, (ns+nc)*N columns)
I_ns = eye(ns);
O_nc = zeros(ns, nc);
A0_dyn = [];
b0_dyn = [];

x0 = [r0; v0; log(m0)];   % initial state in log(m) representation

for j = 1:N
    row = zeros(ns, N*(ns + nc));
    for k = 1:N
        idx_x = (k-1)*(ns+nc) + (1:ns);     % state block positions inside big vector
        idx_u = (k-1)*(ns+nc) + (ns + (1:nc)); % control block positions
        if j == k
            row(:, idx_x) = D(j,k)*I_ns - kt*Ac;
            row(:, idx_u) = -kt*Bc;
        else
            row(:, idx_x) = D(j,k) * I_ns;
            row(:, idx_u) = O_nc;
        end
    end
    A0_dyn = [A0_dyn; row];
    % RHS uses the initial node column (col 1 of D_full) scaled:
    b0_dyn = [b0_dyn; kt*Bc*[g;0]];
end

% Final condition rows (force final state to rf,vf)
A0_fc = zeros(ns-1, N*(ns+nc));
A0_fc(:, (N-1)*ns + (N-1)*nc + (1:ns-1)) = eye(ns-1);
b0_fc = [rf;vf];

A0_ic = zeros(ns, N*ns + N*nc);
A0_ic(:, (1:ns)) = eye(ns); % initial conditions values
b0_ic = [r0;v0; log(m0)];

A0 = [A0_ic;A0_dyn; A0_fc];
b0 = [b0_ic;b0_dyn; b0_fc];

% Pack per-node: [rx; ry; rz; vx; vy; vz; z; ux; uy; uz; s]  -> total per node: ns + nc
% Create final state equality: rx(N)=rf(1), ry(N)=rf(2), rz(N)=rf(3), vx(N)=vf(1), vy(N)=vf(2), vz(N)=vf(3)


%% ------------------ Build linear objective vector c ------------------
% You had c built with quadrature weights;.
c = zeros(N*(ns+nc), 1);
w_l = quad_weights_lobatto(tau);   % quad weights for the full N+1 Radau grid 

for j = 1:N
    % choose index corresponding to the slack variable in node j inside big vector:
    idx_s = (j)*(ns+nc);  
    c(idx_s) = kt * w_l(j);    
end

%% ------------------ CVX Solve ------------------
% Now build the CVX problem. Decision variables are per-node stacked as earlier.
cvx_begin
    cvx_solver sdpt3   % change to mosek if needed
    cvx_precision best

    % Variables: per-node
    variables rx(N) ry(N) rz(N) vx(N) vy(N) vz(N) z(N)
    variables ux(N) uy(N) uz(N) s(N)

    % Stack into a single vector Y to multiply with A0, c, etc.
    Y = [];
    for i=1:N
        Y = [Y; rx(i); ry(i); rz(i); vx(i); vy(i); vz(i); z(i); ux(i); uy(i); uz(i); s(i)];
    end
    % Y = [rx(1); ry(1); rz(1); vx(1); vy(1); vz(1); z(1); ux(1); uy(1); uz(1); s(1)];
    % for i=2:N-1
    %     Y = [Y; rx(i); ry(i); rz(i); vx(i); vy(i); vz(i); z(i); ux(i); uy(i); uz(i); s(i)];
    % end
    % Y = [Y;rx(end); ry(end); rz(end); vx(end); vy(end); vz(end); z(end); ux(end); uy(end); uz(end); s(end)];

    % Objective: minimize c' * Y  
    minimize( c' * Y )
    
    subject to
        % --- Linear dynamics equality and initial and final conditions ---
        A0 * Y == b0;
        %A0_dyn * Y == b0_dyn;

        % --- Thrust cone: ||u|| <= s  (SOC)
        for i = 1:N
            norm([ux(i); uy(i); uz(i)], 2) <= s(i);
        end

        % --- Glideslope convex form: tan(thetaAlt)*sqrt(rx^2+ry^2) <= rz ---
        % enforce rz >= tan(thetaAlt)*sqrt(rx^2 + ry^2)
        for i = 1:N
            norm([rx(i); ry(i)], 2) * tan(thetaAlt) <= rz(i);
            rz(i) >= 0;
        end

        % --- Convexified thrust magnitude bounds  ---
        rho_l = Tl * TmaxTot;
        rho_u = Tu * TmaxTot;
        for i = 1:N
            % (conservative convex envelopes)
            zl = log(m0 - alpha * rho_l * t_nodes(i));
            zu = log(m0 - alpha * rho_u * t_nodes(i));
            % lower approx: rho_l*exp(-zl)*(1 - (z - zl) + 0.5*(z - zl)^2) <= s
            s(i) >= rho_l * exp(-zl) * (1 - (z(i) - zl) + 0.5*(z(i) - zl)^2);

            % upper bound (affine)
            s(i) <= rho_u * exp(-zu) * (1 - (z(i) - zu));
        end
        % for i = 1:N
        %     % lower approx: rho_l*exp(-zl)*(1 - (z - zl) + 0.5*(z - zl)^2) <= s
        %     s(i) >= rho_l * exp(-z(i));
        %     % upper bound (affine)
        %     s(i) <= rho_u * exp(-z(i));
        % end
        
        % for i= 1:N
        %     z_bar = [z(i);s(i)];
        % 
        %     %zl = log(m0-alpha*rho_l*t_nodes(i+1)); zu = log(m0-alpha*rho_u*t_nodes(i+1));
        %     zl = log(m0-alpha*rho_l*t_nodes(i)); zu = log(m0-alpha*rho_u*t_nodes(i));
        %     A_rho = [rho_l*exp(-zl/sqrt(2)), 0];
        %     %A_rho = [0.5*rho_l*exp(-zl), 0];
        %     b_rho = -[rho_l*exp(-zl)+zl, 1];
        %     %b_rho = -[rho_l*exp(-zl)*(1+zl), 1];
        %     c_rho = double(rho_l*exp(-zl)*(1+zl+0.5*zl^2));
        % 
        %     A3_c = [0.5*b_rho; A_rho];
        %     b3_c = [0.5*c_rho + 0.5; 0];
        %     c3_c = -b_rho*0.5;
        %     d3_c = 0.5-c_rho*0.5;
        % 
        %     A4_c = [rho_u*exp(-zu), 1];
        %     c4_c = rho_u*exp(-zu)*(1+zu);
        % 
        %     norm(A3_c*z_bar + b3_c,2) <= (c3_c*z_bar + d3_c);
        %     %A3_c*z_bar + b3_c <= (c3_c*z_bar + d3_c)
        %     A4_c*z_bar <= c4_c*z_bar;
        %     %A4_c*z_bar <= c4_c;
        %     %norm(A4_c*z_bar,2) <= c4_c*z_bar;
        %     %norm(A4_c*z_bar,2) <= c4_c;
        % 
        % end

        % --- box and sign constraints ---
        % s >= 0;
        % z <= log(m0);   % log mass cannot exceed initial (optional)
        % z >= log(m_dry);
        % keep velocities, positions free (boundaries handled below)

        % % --- Boundary conditions: initial log-mass and final r/v ---
        % % only enforce final state equality at collocation node N (which corresponds to last interior node).
        % rx(end) == rf(1);
        % ry(end) == rf(2);
        % rz(end) == rf(3);
        % 
        % vx(end) == vf(1);
        % vy(end) == vf(2);
        % vz(end) == vf(3);
        % 
        % rx(1) == r0(1);
        % ry(1) == r0(2);
        % rz(1) == r0(3);
        % 
        % vx(1) == v0(1);
        % vy(1) == v0(2);
        % vz(1) == v0(3);
        % 
        % z(1) == log(m0);

cvx_end

%% ------------------ Unpack & postprocess ------------------
% Build full state arrays including the initial state r0,v0,m0 for plotting (as you did)
r_sol = [rx(:)'; ry(:)'; rz(:)']';
v_sol = [vx(:)'; vy(:)'; vz(:)']';
z_sol = z(:);
m_sol = [exp(z_sol)];
u_sol = [ux(:), uy(:), uz(:)];
sigma_sol = s(:);
% Throttle_sol=zeros(size(sigma_sol));
% for j=1:N
%     Throttle_sol(j) = norm([ux(j) uy(j) uz(j)])/6;
% end 
Tnorm_sol = m_sol.*sigma_sol;%m_sol.*Throttle_sol;

fprintf('CVX status: %s\n', cvx_status);
% if exist('cvx_optval','var')
%     fprintf('Objective value (c''*Y) = %.6g\n', cvx_optval);
% end

fprintf('Mass consumed: %.3f kg\n', abs(m_sol(end)-m0));

%% ------------------ Plots ------------------
t_plot = t_nodes(:);

% Position plot
figure(1);

subplot(3,1,1);
plot(t_plot, r_sol(:,1), 'r-','LineWidth',1.6); hold on;
plot(t_plot, r_sol(:,2), 'g--','LineWidth',1.6);
plot(t_plot, r_sol(:,3), 'b-.','LineWidth',1.6);
xlabel('t (s)'); ylabel('r (m)'); legend('x','y','z'); grid on; title('Position');
set(gca, 'FontSize', 18);

% Velocity plot
subplot(3,1,2);
plot(t_plot, v_sol(:,1), 'r-','LineWidth',1.6); hold on;
plot(t_plot, v_sol(:,2), 'g--','LineWidth',1.6);
plot(t_plot, v_sol(:,3), 'b-.','LineWidth',1.6);
xlabel('t (s)'); ylabel('v (m/s)'); legend('vx','vy','vz'); grid on; title('Velocity');
set(gca, 'FontSize', 18);

% Mass & thrust
subplot(3,1,3);
yyaxis left; plot(t_plot, m_sol,'-','LineWidth',1.6); ylabel('Mass [kg]');
yyaxis right; plot(t_nodes, Tnorm_sol/TmaxTot,'-','LineWidth',1.4); ylabel('||T||/m ');
xlabel('t (s)'); grid on; title('Mass & Throttle');
set(gca, 'FontSize', 18);

saveas(gcf, 'Results/cvx-opt_mars_descent_results-lobatto.png');


% 3D Trajectory with thrust vectors
figure(2); clf;
plot3(r_sol(:,1), r_sol(:,2), r_sol(:,3), 'b-','LineWidth',1.6); hold on;
plot3(r0(1), r0(2), r0(3), 'go', 'MarkerSize',8, 'MarkerFaceColor','g');
plot3(rf(1), rf(2), rf(3), 'rx', 'MarkerSize',10, 'LineWidth',1.6);

% Scale factor for thrust arrows (adjust for visibility)
scale = 0.3;  

% Plot thrust vectors at every k-th node to avoid clutter
k_skip = 1; 
quiver3(r_sol(1:k_skip:end,1), r_sol(1:k_skip:end,2), r_sol(1:k_skip:end,3), ...
        ux(1:k_skip:end), uy(1:k_skip:end), uz(1:k_skip:end), ...
        scale, 'r', 'LineWidth',1.2, 'MaxHeadSize',2);

grid on; axis equal; view(3);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('3D Trajectory');
legend('Trajectory','Start','Target','Thrust direction');

set(gca, 'FontSize', 18);
saveas(gcf, 'Results/cvx-opt_mars_descent_trajectory-lobatto.png');



% Save
save('mars_traj_cvx_sol.mat', 'rx','ry','rz','vx','vy','vz','z','ux','uy','uz','s','m_sol','Tnorm_sol','t_nodes');

fprintf('Saved solution to mars_traj_cvx_sol-lobatto.mat\n');


%% Animated Descent Plot
figure(3); clf;
hold on; grid on;
axis equal;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('Mars Descent Animation');
set(gca, 'FontSize', 18);
view(3);

% Plot terrain reference
terrain = fill3([-1000, 1000, 1000, -1000]*3, ...
                [-1000, -1000, 1000, 1000]*3, ...
                [0, 0, 0, 0], [0.5 0.3 0.1], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName','Ground');
            
traj_line = plot3(NaN, NaN, NaN, 'b-', 'LineWidth', 2, 'DisplayName','Trajectory');
lander = plot3(NaN, NaN, NaN, 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r', 'DisplayName','Lander');

time_txt = text(r0(1), r0(2), r0(3)+100, '', 'FontSize', 12, 'FontWeight', 'bold');

legend;

% --- Animation and Video Writer ---
v = VideoWriter('Results/mars_descent_cvx_lobattto.mp4', 'MPEG-4');
v.Quality = 100;
v.FrameRate = 10;
open(v);
writeVideo(v, getframe(gcf));

% Animation loop
for k = 1:length(t_nodes)
    traj_line.XData = r_sol(1:k,1);
    traj_line.YData = r_sol(1:k,2);
    traj_line.ZData = r_sol(1:k,3);
    
    lander.XData = r_sol(k,1);
    lander.YData = r_sol(k,2);
    lander.ZData = r_sol(k,3);
    
    time_txt.Position = [r_sol(k,1), r_sol(k,2), r_sol(k,3)+100];
    time_txt.String = sprintf('t = %.1f s', t_nodes(k));
    
    drawnow;
    writeVideo(v, getframe(gcf));
end
writeVideo(v, getframe(gcf));
close(v);
disp('Animation saved to trajectory_animation.mp4');

traj_line.XData = r_sol(1:end,1);
traj_line.YData = r_sol(1:end,2);
traj_line.ZData = r_sol(1:end,3);

lander.XData = r_sol(end,1);
lander.YData = r_sol(end,2);
lander.ZData = r_sol(end,3);

time_txt.Position = [r_sol(end,1), r_sol(end,2), r_sol(end,3)+100];
time_txt.String = sprintf('t = %.1f s', t_nodes(end));

pause(0.2);
drawnow;


%% Comparison plots with non conex 
data = load('mars_traj_correct.mat');
x_0 = data.xsol;   % Use as initial guess
N0 = 30;
tau0 = lobatto_nodes(N0);
t_0 = kt*tau0 + (tf+t0)/2; 
r_0 = reshape(x_0(1:3*N0), N0, 3);
v_0 = reshape(x_0(3*N0+1:6*N0), N0, 3);
m_0 = x_0(6*N0+1 : 7*N0);
T_0 = reshape(x_0(7*N0+1:10*N0), N0, 3);


% Time Array
tau = lobatto_nodes(N);
t = kt*tau + (tf+t0)/2; 

xsol0 = lagrange_interp_matrix_barycentric(t_0, r_0(:,1), t);
ysol0 = lagrange_interp_matrix_barycentric(t_0, r_0(:,2), t);
zsol0 = lagrange_interp_matrix_barycentric(t_0, r_0(:,3), t);
usol0 = lagrange_interp_matrix_barycentric(t_0, v_0(:,1), t);
vsol0 = lagrange_interp_matrix_barycentric(t_0, v_0(:,2), t);
wsol0 = lagrange_interp_matrix_barycentric(t_0, v_0(:,3), t);
msol0 = lagrange_interp_matrix_barycentric(t_0, m_0, t);
T1sol0 = lagrange_interp_matrix_barycentric(t_0, T_0(:,1), t);
T2sol0 = lagrange_interp_matrix_barycentric(t_0, T_0(:,2), t);
T3sol0 = lagrange_interp_matrix_barycentric(t_0, T_0(:,3), t);
X0 = [xsol0;ysol0;zsol0;usol0;vsol0;wsol0;msol0;T1sol0;T2sol0;T3sol0];

figure(100); clf;
plot(t_plot, m_sol, 'b-', 'LineWidth', 1.6, 'DisplayName','CVX'); hold on;
plot(t, msol0, 'r--', 'LineWidth', 1.6, 'DisplayName','Non-CVX');

xlabel('Time [s]');
ylabel('Mass [kg]');
title('Mass Consumption Comparison');
grid on;
legend('Location','best');
set(gca, 'FontSize', 18);

saveas(gcf, 'Results/opt_mars_descent_mass_compareL.png');



figure(101); clf;
plot(m_sol, rz,  'b-', 'LineWidth', 1.6, 'DisplayName','CVX'); hold on;
plot(msol0, zsol0, 'r--', 'LineWidth', 1.6, 'DisplayName','Non-CVX');

ylabel('Height [m]');
xlabel('Mass [kg]');
title('Mass Consumption Comparison with Height');
grid on;
legend('Location','best');
set(gca, 'FontSize', 18);

saveas(gcf, 'Results/opt_mars_descent_mass_compareL_z.png');