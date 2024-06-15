% Parameters
num_links = 5; % Number of links
link_length = 1; % Length of each link
Ts = 0.1; % Sample time
T = 10; % Total simulation time
time = 0:Ts:T; % Time vector

% Initial joint angles (all set to zero)
theta = zeros(num_links, length(time));

% Joint velocities (random for demonstration)
joint_velocities = randn(num_links, length(time)) * 0.1;

% Compute the positions of each joint over time
positions = zeros(2, num_links+1, length(time));

for t = 2:length(time)
    % Update joint angles
    theta(:, t) = theta(:, t-1) + joint_velocities(:, t) * Ts;
    
    % Compute positions using forward kinematics
    for i = 1:num_links
        positions(:, i+1, t) = positions(:, i, t) + ...
            [cos(sum(theta(1:i, t))); sin(sum(theta(1:i, t)))] * link_length;
    end
end

% Plot the results
figure;
for t = 1:length(time)
    plot(positions(1, :, t), positions(2, :, t), 'o-', 'LineWidth', 2);
    axis equal;
    axis([-num_links, num_links, -num_links, num_links]);
    title(sprintf('Time: %.1f s', time(t)));
    xlabel('X Position');
    ylabel('Y Position');
    drawnow;
    pause(Ts);
end
