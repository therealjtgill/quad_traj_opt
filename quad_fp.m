clear all; clc;
%% https://www.mathworks.com/help/matlab/creating_plots/animation-techniques-1.html
%% Unpack named states from state vector
function [pos_W, vel_W, rpy_rad, ang_vel_B] = unpack_state(x)
    pos_W = x(1:3);
    vel_W = x(4:6);
    rpy_rad = x(7:9);
    ang_vel_B = x(10:12);
end

%% Pack named states into one state vector
function x = pack_states(pos_W, vel_W, rpy_rad, ang_vel_B)
    x = zeros(12, 1);
    x(1:3) = pos_W;
    x(4:6) = vel_W;
    x(7:9) = rpy_rad;
    x(10:12) = ang_vel_B;
end

%% Unpack state and input from decision vector
function [X, U] = unpack(z, nx, nu)
    % Horizon length
    N = (length(z) - nx) / (nx + nu); % Exclude final state
    X = reshape(z(1: nx*(N + 1)), nx, N + 1);
    U = reshape(z(nx * (N + 1) + 1:end), nu, N);
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

% % Moment of inertia tensor, kind of a guess
% J = diag([0.5, 0.5, 1.0]);

%% Quadrotor dynamics
function x_dot = quad_dynamics(x, u)
    hover_alpha = 0.3;
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
    motor_speeds = u * max_motor_speed_radps;
    R_B_to_W = rot_B_to_W(rpy_rad(1), rpy_rad(2), rpy_rad(3));

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

    % Angular acceleration in body frame
    ang_vel_B_dot = inv(J) * (cross(-ang_vel_B, J * ang_vel_B) + torque_B - c_ang_drag * ang_vel_B * norm(ang_vel_B));

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

%% Constraints
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
    last_state_dynamics = quad_dynamics(X(:, end - 1), U(:, end));
    ceq = [ceq; last_state_dynamics(4:6)];

    % % Trapezoidal implicit integration
    % for k = 1:N
    %     f_k = quad_dynamics(X(:, k), U(:, k));
    %     f_kp1 = quad_dynamics(X(:, k + 1), U(:, k));
    % 
    %     ceq = [ceq; X(:, k + 1) - X(:, k) - 0.5 * dt * (f_k + f_kp1)];
    % end

    for k = 1:N
        k1 = quad_dynamics(X(:, k), U(:, k));
        k2 = quad_dynamics(X(:, k) + k1 * dt / 2, U(:, k));
        k3 = quad_dynamics(X(:, k) + k2 * dt / 2, U(:, k));
        k4 = quad_dynamics(X(:, k) + k3 * dt, U(:, k));

        ceq = [ceq; X(:, k + 1) - X(:, k) - (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)];
    end

    % Individual motor commands must be between 0 and 1 (motors can't spin
    % backwards)
    cineq = [-U; U - 1];
end


%% Optimization parameters
N = 45;
T = 3;
dt = T / N;

nx = 12;
nu = 4;

% Boundary conditions
x0 = zeros(12, 1);
xT = zeros(12, 1);
xT(1) = 0;
xT(2) = 0;
xT(3) = 0;
xT(8) = pi*2;
xT(9) = 0;

z0 = zeros(nx * (N + 1) + nu * N, 1);

% Initial guess: linear interpolation
for k = 0:N
    z0(nx*k + 1: nx*(k + 1)) = x0 + (k/N)*(xT-x0);
end

% Position initialization just for the 360 flip
for k = 0:N
    if k < 2N
        z0(nx*k + 1: nx*k + 3) = [0 0 10 * 2 * k / N]';
    else
        z0(nx*k + 1: nx*k + 3) = [0 0 10 * 2 * k / N - 10 * (k - N) / (2 * N)]';
    end
end

% Set control input to half throttle for all motors
z0(nx*(k + 1) + 1:end) = 0.5;

%% Solve NLP and make plots

% Optimizer settings
options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','iter', ...
    'MaxFunctionEvaluations',3e5);
% options = optimoptions('fmincon', ...
%     'Algorithm','interior-point', ...
%     'EnableFeasibilityMode', true, ...
%     'Display','iter', ...
%     'SubproblemAlgorithm', 'cg', ...
%     'MaxFunctionEvaluations',1e5);

tic

[z, cost_val] = fmincon(@(z) cost(z, nu, nu, T), z0, [], [], [], [], [], [], @(z) constraints(z, nx, nu, T, x0, xT), options);
t = toc;

% Extract state and input trajectories
[X,U] = unpack(z,nx,nu);
% Plotting
plot_results(X,U,T)
% subplot(4,1,1)
title(['Collocation, cost = ' num2str(cost_val) ', time = ' num2str(t) '[s]'],FontSize=18)

