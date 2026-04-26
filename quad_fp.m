clear all; clc;

%% Unpack named states from state vector
function [pos_W, vel_W, rpy_rad, ang_vel_B] = unpack_state(x)
    pos_W = x(1:3);
    vel_W = x(4:6);
    rpy_rad = x(7:9);
    ang_vel_B = x(10:12);
end

%% Pack named states into one state vector
function x = pack_states(pos_W, vel_W, rpy_rad, ang_vel_B)
    % disp(size(ang_vel_B))
    x = zeros(12, 1);
    x(1:3) = pos_W;
    x(4:6) = vel_W;
    x(7:9) = rpy_rad;
    x(10:12) = ang_vel_B;
end

%% Unpack state and input from decision vector
function [X, U] = unpack(z, nx, nu)
    % Computing horizon length
    N = (length(z) - nx) / (nx + nu); % Substracting final state cost first

    X = reshape(z(1: nx*(N + 1)), nx, N + 1);
    % U = z(nx * (N + 1) + 1:end);
    U = reshape(z(nx * (N + 1) + 1:end), nu, N);
end

%% Rotate from body frame to world frame (FRD to NED matrix from RPY in radians)
function R_B_to_W = rot_B_to_W(roll, pitch, yaw)
    sphi = sin(roll);
    cphi = cos(roll);
    stheta = sin(pitch);
    ctheta = cos(pitch);
    spsi = sin(yaw);
    cpsi = cos(yaw);

    R_B_to_W = zeros(3,3);
    R_B_to_W(1, 1) = ctheta * cpsi;
    R_B_to_W(1, 2) = sphi * stheta * cpsi - cphi * spsi;
    R_B_to_W(1, 3) = cphi * stheta * cpsi + sphi * spsi;
    
    R_B_to_W(2, 1) = ctheta * spsi;
    R_B_to_W(2, 2) = sphi * stheta * spsi + cphi * cpsi;
    R_B_to_W(2, 3) = cphi * stheta * spsi - sphi * cpsi;
    
    R_B_to_W(3, 1) = -stheta;
    R_B_to_W(3, 2) = sphi * ctheta;
    R_B_to_W(3, 3) = cphi * ctheta;
end

%% Quadrotor parameters
% hover_alpha = 0.5;
% max_motor_speed_radps = 3000.0;
% min_motor_speed_radps = 0.0;
% 
% mass_kg = 0.25;
% k_thrust = mass_kg * 9.81 / (4.0 * (hover_alpha * max_motor_speed_radps)^2);
% c_lin_drag = 0.2;
% k_torque = 2.163 / 100000;
% c_ang_drag = 0.1;
% 
% % Each column is the location of the end of an arm in body frame
% arms_B = [
%     0.125  0.125 -0.125 -0.125;
%     0.125 -0.125 -0.125  0.125;
%     0.0      0.0    0.0    0.0;
% ];
% 
% % Moment of inertia tensor, kind of a guess
% J = diag([0.5, 0.5, 1.0]);

%% Quadrotor dynamics
function x_dot = quad_dynamics(x, u)
    % global max_motor_speed_radps;
    % global min_motor_speed_radps;
    % global mass_kg;
    % global k_thrust;
    % global c_lin_drag;
    % global k_torque;
    % global c_ang_drag;
    % global arms_B;
    % global J;

    hover_alpha = 0.5;
    max_motor_speed_radps = 3000.0;
    min_motor_speed_radps = 0.0;
    
    mass_kg = 0.25;
    k_thrust = mass_kg * 9.81 / (4.0 * (hover_alpha * max_motor_speed_radps)^2);
    c_lin_drag = 0.2;
    k_torque = 2.163 / 100000;
    c_ang_drag = 0.1;
    
    % Each column is the location of the end of an arm in body frame
    arms_B = [
        0.125  0.125 -0.125 -0.125;
        0.125 -0.125 -0.125  0.125;
        0.0      0.0    0.0    0.0;
    ];
    
    % Moment of inertia tensor, kind of a guess
    J = diag([0.5, 0.5, 1.0]);

    [pos_W, vel_W, rpy_rad, ang_vel_B] = unpack_state(x);
    % disp(size(u))
    motor_speeds = u * max_motor_speed_radps;
    R_B_to_W = rot_B_to_W(rpy_rad(1), rpy_rad(2), rpy_rad(3));

    % x_dot = zeros(12, 1);
    % [pos_W_dot, vel_W_dot, rpy_rad_dot, ang_vel_B_dot] = unpack_state(x_dot);
    pos_W_dot = zeros(3, 1);
    vel_W_dot = zeros(3, 1);
    rpy_rad_dot = zeros(3, 1);
    ang_vel_B_dot = zeros(3, 1);

    % Gravity and quadratic velocity drag
    vel_W_dot = (1.0 / mass_kg) * ([0; 0; -9.81] - c_lin_drag * vel_W * norm(vel_W));

    % Motor forces in world frame
    motor_force = zeros(3, 1);
    for i = 1:4
        motor_force = motor_force + R_B_to_W * ([0 0 k_thrust * motor_speeds(i)^2]');
    end
    vel_W_dot = vel_W_dot + motor_force / mass_kg;

    pos_W_dot = vel_W;

    torque_B = zeros(3, 1);

    % Body frame torque from spinning motors
    for i = 1:4
        torque_B = torque_B + (2 * (mod(i, 2) == 0) - 1) * k_torque * motor_speeds(i) ^ 2 * [0 0 1]';
    end

    % Body frame torque from applying forces on levers about the CM
    for i = 1:4
        thrust_B_i = (1 / mass_kg) * [0; 0; k_thrust * motor_speeds(i)^2];
        torque_B = torque_B + cross(arms_B(:, i), thrust_B_i);
    end

    ang_vel_B_dot = inv(J) * (cross(-ang_vel_B, J * ang_vel_B) + torque_B - c_ang_drag * ang_vel_B * norm(ang_vel_B));
    % disp("main cross product")
    % disp(size(cross(-ang_vel_B, J * ang_vel_B)))
    % disp("torque_b")
    % disp(size(torque_B))

    rpy_rad_dot(1) = ang_vel_B(1) + sin(rpy_rad(1)) * tan(rpy_rad(2)) * ang_vel_B(2) + cos(rpy_rad(1)) * tan(rpy_rad(2)) * ang_vel_B(3);
    rpy_rad_dot(2) = cos(rpy_rad(1)) * ang_vel_B(2) - sin(rpy_rad(1)) * ang_vel_B(3);
    rpy_rad_dot(3) = (sin(rpy_rad(1) / cos(rpy_rad(2)))) * ang_vel_B(2) + (cos(rpy_rad(1) / cos(rpy_rad(2)))) * ang_vel_B(3);

    x_dot = pack_states(pos_W_dot, vel_W_dot, rpy_rad_dot, ang_vel_B_dot);
end

%% Cost function
function J = cost(z, nx, nu, T)
    % z:    decision variables
    % nx:   size of state variables
    % nu:   size of input variable
    % T:    final time

    [~, U] = unpack(z, nx, nu);
    U_flat = reshape(U, 1, []);
    N = length(U);
    dt = T / N;

    J = dt * sum(U_flat.^2);
end

%% Constaints
function [cineq, ceq] = constraints(z, nx, nu, T, x0, xT)
    % z:    decision variables
    % nx:   size of state variables
    % nu:   size of input variable
    % T:    final time
    % x0:   initial state
    % xT:   desired final state

    [X, U] = unpack(z, nx, nu);

    N = length(U);
    dt = T/N;

    ceq = [];
    ceq = [ceq; X(:, 1) - x0; X(:, end) - xT];

    for k = 1:N
        f_k = quad_dynamics(X(:, k), U(:, k));
        f_kp1 = quad_dynamics(X(:, k + 1), U(:, k));

        ceq = [ceq; X(:, k + 1) - X(:, k) - 0.5 * dt * (f_k + f_kp1)];
    end

    % Individual motor commands must be between 0 and 1
    cineq = [-U; U - 1.0];
end

%% Plotting function
function plot_results(X,U,T)
    N = length(U);
    tX = linspace(0,T,N+1);
    tU = linspace(0,T,N);
    
    figure
    subplot(4,1,1)
    plot(tX,X(1,:),'LineWidth',3); grid on
    ylabel('$x$','Interpreter','latex')
    
    subplot(4,1,2)
    plot(tX,X(2,:),'LineWidth',3); grid on
    ylabel('$y$','Interpreter','latex')
    
    subplot(4,1,2)
    plot(tX,X(3,:),'LineWidth',3); grid on
    ylabel('$z$','Interpreter','latex')
    
    subplot(4,1,4)
    stairs(tU,U(1, :),'LineWidth',3); grid on;
    ylabel('$u$','Interpreter','latex')
    xlabel('time')
    
    fontsize(18,"points")
end

%% Optimization parameters
N = 100;
T = 3;
dt = T / N;

nx = 12;
nu = 4;

% Boundary conditions
x0 = zeros(12, 1);
xT = zeros(12, 1);
xT(1) = 1;
xT(2) = 1;
xT(3) = 1;

z0 = zeros(nx * (N + 1) + nu * N, 1);

% Initial guess: linear interpolation
for k = 0:N
    z0(nx*k+1:nx*(k+1)) = x0 + (k/N)*(xT-x0);
end

%% Solve NLP

% Optimizer settings
options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','iter', ...
    'MaxFunctionEvaluations',5e5);

tic

[z, cost_val] = fmincon(@(z) cost(z, nu, nu, T), z0, [], [], [], [], [], [], @(z) constraints(z, nx, nu, T, x0, xT), options);
t = toc;

% Extracting state and input trajectories
[X,U] = unpack(z,nx,nu);
% Plotting
plot_results(X,U,T)
subplot(3,1,1)
title(['Collocation, cost = ' num2str(cost_coll) ', time = ' num2str(t) '[s]'],FontSize=18)

