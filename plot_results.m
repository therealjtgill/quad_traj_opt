%% Plotting function
function plot_results(X,U,T)
    N = length(U);
    tX = linspace(0,T,N+1);
    tU = linspace(0,T,N);
    
    % Motor commands
    figure(1)
    subplot(4,1,1)
    stairs(tU,U(1, :),'LineWidth',3); grid on;
    ylabel('$u$','Interpreter','latex')
    xlabel('time')

    subplot(4,1,2)
    stairs(tU,U(2, :),'LineWidth',3); grid on;
    ylabel('$u$','Interpreter','latex')
    xlabel('time')

    subplot(4,1,3)
    stairs(tU,U(3, :),'LineWidth',3); grid on;
    ylabel('$u$','Interpreter','latex')
    xlabel('time')

    subplot(4,1,4)
    stairs(tU,U(4, :),'LineWidth',3); grid on;
    ylabel('$u$','Interpreter','latex')
    xlabel('time')

    % Each column is the location of the end of an arm in body frame
    arms_B = [
        0.125  0.125 -0.125 -0.125;
        0.125 -0.125 -0.125  0.125;
        0.0      0.0    0.0    0.0;
    ];

    % Sampled drone locations
    figure(2)
    % plot3(X(1,:), X(2,:), X(3,:), 'LineWidth', 3); grid on; hold on;

    positions = X(1:3, :);
    rpy_rads = X(7:9, :);
    m1_positions = [];
    m2_positions = [];
    m3_positions = [];
    m4_positions = [];
    for k = 1:1:(N + 1)
        rpy_rad = rpy_rads(:, k);
        R_B_to_W = rot_B_to_W(rpy_rad(1), rpy_rad(2), rpy_rad(3));
        m1_positions = [m1_positions, R_B_to_W * arms_B(:, 1) + positions(:, k)];
        m2_positions = [m2_positions, R_B_to_W * arms_B(:, 2) + positions(:, k)];
        m3_positions = [m3_positions, R_B_to_W * arms_B(:, 3) + positions(:, k)];
        m4_positions = [m4_positions, R_B_to_W * arms_B(:, 4) + positions(:, k)];
    end

    plot3(m1_positions(1, :), m1_positions(2, :), m1_positions(3, :), '--'); hold on;
    plot3(m2_positions(1, :), m2_positions(2, :), m2_positions(3, :), '--'); hold on;
    plot3(m3_positions(1, :), m3_positions(2, :), m3_positions(3, :), '--'); hold on;
    plot3(m4_positions(1, :), m4_positions(2, :), m4_positions(3, :), '--'); hold on;

    scatter3(m1_positions(1, :), m1_positions(2, :), m1_positions(3, :)); hold on;
    scatter3(m2_positions(1, :), m2_positions(2, :), m2_positions(3, :)); hold on;
    scatter3(m3_positions(1, :), m3_positions(2, :), m3_positions(3, :)); hold on;
    scatter3(m4_positions(1, :), m4_positions(2, :), m4_positions(3, :)); 

    for k = 1:1:(N + 1)
        rpy_rad = rpy_rads(:, k);
        R_B_to_W = rot_B_to_W(rpy_rad(1), rpy_rad(2), rpy_rad(3));
        m1_position = R_B_to_W * arms_B(:, 1) + positions(:, k);
        arm_line = [m1_position, positions(:, k)];
        plot3(arm_line(1, :), arm_line(2, :), arm_line(3, :), 'black');
        m2_position = R_B_to_W * arms_B(:, 2) + positions(:, k);
        arm_line = [m2_position, positions(:, k)];
        plot3(arm_line(1, :), arm_line(2, :), arm_line(3, :), 'black');
        m3_position = R_B_to_W * arms_B(:, 3) + positions(:, k);
        arm_line = [m3_position, positions(:, k)];
        plot3(arm_line(1, :), arm_line(2, :), arm_line(3, :), 'black');
        m4_position = R_B_to_W * arms_B(:, 4) + positions(:, k);
        arm_line = [m4_position, positions(:, k)];
        plot3(arm_line(1, :), arm_line(2, :), arm_line(3, :), 'black');
    end

    % Plot everything inside of a cube that's guaranteed to enclose the
    % entire trajectory
    xs = X(1, :);
    ys = X(2, :);
    zs = X(3, :);

    dx = max(xs) - min(xs);
    dy = max(ys) - min(ys);
    dz = max(zs) - min(zs);

    box_dim = max([dx dy dz]) * 1.05;

    box_center = [(max(xs) + min(xs)) / 2; (max(ys) + min(ys)) / 2; (max(zs) + min(zs)) / 2];

    xlim([box_center(1) - box_dim/2, box_center(1) + box_dim/2]);
    ylim([box_center(2) - box_dim/2, box_center(2) + box_dim/2]);
    zlim([box_center(3) - box_dim/2, box_center(3) + box_dim/2]);

    ylabel('drone position over time')

    fontsize(18,"points")
end
