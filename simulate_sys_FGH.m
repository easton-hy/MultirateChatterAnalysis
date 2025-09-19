clear all
close all
clc
% Add existing paths
addpath('data','src');
% Use existing parameter setup
make_Global
make_SDM
TOOL = 1; % Define tool type (1, 2, 3, 4, or 5 - see make_Tool.m)
o = 2800 * 2*pi/60; % spindle speed (rad/s)
make_Tool
make_Work
% Extract M, C, K from the tool parameters (for TOOL=1)
% From make_Tool.m, TOOL=1 creates scalar values
M_scalar = mass; % 0.4 kg
K_scalar = K; % Stiffness from make_Tool.m
B_scalar = B; % Damping from make_Tool.m
% Create 2x2 matrices as in the paper's formulation
M_matrix = M_scalar * eye(2); % Mass matrix
C_matrix = B_scalar * eye(2); % Damping matrix
K_matrix = K_scalar * eye(2); % Stiffness matrix
% Simulation parameters for Eq. 7
tau = 2*pi/o/FLUTES; % time delay
dt = tau/2500; % time step
T_sim = 0.03; % simulation time (s)
t = 0:dt:T_sim;
N = length(t);
% Calculate Kf matrix using existing function
Work.Kf = calc_Kf_CS(Tool,SDM,Work);
% Initialize state vectors for Eq. 7
% x = [q; q_dot] = [x_pos; y_pos; x_vel; y_vel]
x = zeros(4, N);
x(:,1) = [0; 0; 0; 0]; % initial condition

% Initialize acceleration storage
y_acc = zeros(1, N);

% Feed input
uf = 1e-5; % feed per tooth (m)
u = uf * ones(N,1);
% Axial depth of cut - will be modified during simulation
w_initial = 1.5e-3; % 1.5 mm (this is ap in the paper)
fprintf('Simulating Equation 7...\n');
% Main simulation loop
for i = 2:N
    % Current time and spindle angle
    t_curr = t(i);
    theta = mod(o * t_curr, 2*pi);
    % Set axial depth of cut to 0 after t > 0.03s
    if t_curr > 0.03
        ap = 0; % No cutting after 0.03s
    else
        ap = w_initial; % Normal cutting before 0.03s
    end
    % Determine which Kf to use based on current angle
    k_idx = floor(theta / (2*pi) * SDM.k) + 1;
    if k_idx > SDM.k
        k_idx = SDM.k;
    end
    % Current Kf matrix (2x2)
    Kf_curr = Work.Kf(:,:,k_idx);
    fprintf('Kf matrix values:\n');
    disp(Kf_curr);
    fprintf('Max Kf: %.2e\n', max(Kf_curr(:)));
    % Form A1(t), A2(t), B(t) matrices EXACTLY as in Eq. 7 (equations 10-12)
    % A1(t) = [0, I; -M^(-1)(K + ap*Kf(t)), -M^(-1)*C]
    A1 = [zeros(2,2), eye(2);
          -inv(M_matrix)*(K_matrix + ap*Kf_curr), -inv(M_matrix)*C_matrix];
    % A2(t) = [0, 0; ap*M^(-1)*Kf(t), 0]
    A2 = [zeros(2,2), zeros(2,2);
          ap*inv(M_matrix)*Kf_curr, zeros(2,2)];
    % B(t) = [0; ap*M^(-1)*Kf(t)] for feed in x-direction
    B_eq7 = [zeros(2,1);
             ap*inv(M_matrix)*Kf_curr*[1; 0]];
    % Calculate delayed state x(t-tau)
    delay_steps = round(tau/dt);
    if i > delay_steps
        x_delayed = x(:, i-delay_steps);
    else
        x_delayed = zeros(4, 1);
    end
    % Equation 7: x'(t) = A1(t)x(t) + A2(t)x(t-tau) + B(t)u(t)
    x_dot = A1 * x(:,i-1) + A2 * x_delayed + B_eq7 * u(i);
    x(:,i) = x(:,i-1) + dt * x_dot;
    
    % Extract y-acceleration directly from x_dot (4th component)
    y_acc(i) = x_dot(4);  % y-acceleration is d²y/dt²
end
fprintf('Simulation completed.\n');

% Extract positions and velocities
q_pos = x(1:2,:); % [x_pos; y_pos]
q_vel = x(3:4,:); % [x_vel; y_vel]

% Plot acceleration over time
figure;
plot(t, y_acc, 'r-', 'LineWidth', 1.5);  % Convert to milliseconds
xlabel('Time (ms)');
ylabel('Y Acceleration (m/s²)');
title(sprintf('Y-Acceleration vs Time (ω = %.0f rpm) - Following Eq. 7', o*60/(2*pi)));
grid on;
% Add vertical line at t = 0.03s to show when cutting stops
hold on;
xline(0.2, 'k--', 'Cutting stops', 'LineWidth', 1);  % 30 ms = 0.03 s
% Display system parameters for verification
fprintf('\nSystem Parameters:\n');
fprintf('Mass: %.3f kg\n', M_scalar);
fprintf('Stiffness: %.2e N/m\n', K_scalar);
fprintf('Damping: %.1f N⋅s/m\n', B_scalar);
fprintf('Natural frequency: %.1f Hz\n', sqrt(K_scalar/M_scalar)/(2*pi));