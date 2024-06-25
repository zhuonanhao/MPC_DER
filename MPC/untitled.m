clear; clc; close all;

%% System Parameters
l1 = 1.0;  % Length of link 1
l2 = 1.0;  % Length of link 2
m1 = 1.0;  % Mass of link 1
m2 = 1.0;  % Mass of link 2
I1 = 0.1;  % Inertia of link 1
I2 = 0.1;  % Inertia of link 2
g = 9.81;  % Acceleration due to gravity

%% Dynamics Functions
% Mass (Inertia) Matrix
M = @(theta) [I1 + I2 + m1*(l1/2)^2 + m2*(l1^2 + (l2/2)^2 + 2*l1*(l2/2)*cos(theta(2))), I2 + m2*((l2/2)^2 + l1*(l2/2)*cos(theta(2)));
              I2 + m2*((l2/2)^2 + l1*(l2/2)*cos(theta(2))), I2 + m2*(l2/2)^2];

% Coriolis and Centrifugal Matrix
C = @(theta, dtheta) [-m2*l1*(l2/2)*sin(theta(2))*dtheta(2), -m2*l1*(l2/2)*sin(theta(2))*(dtheta(1) + dtheta(2));
                      m2*l1*(l2/2)*sin(theta(2))*dtheta(1), 0];

% Gravity Vector
G = @(theta) [m1*g*(l1/2)*cos(theta(1)) + m2*g*(l1*cos(theta(1)) + (l2/2)*cos(theta(1) + theta(2)));
              m2*g*(l2/2)*cos(theta(1) + theta(2))];

% Dynamics
dynamics = @(t, x, u) [x(3:4);
                       M(x(1:2)) \ (u - C(x(1:2), x(3:4))*x(3:4) - G(x(1:2)))];

%% Initial Conditions and Simulation Parameters
theta1_0 = 0;  % Initial angle of joint 1
theta2_0 = 0;  % Initial angle of joint 2
dtheta1_0 = 0;    % Initial angular velocity of joint 1
dtheta2_0 = 0;    % Initial angular velocity of joint 2

x0 = [theta1_0; theta2_0; dtheta1_0; dtheta2_0];  % Initial state

T = 10;  % Simulation time
dt = 0.01;  % Time step
time = 0:dt:T;

%% Control Input (for simplicity, we use zero torques here)
u = @(t) [0; 0.5];

%% Simulate the System
[t, x] = ode45(@(t, x) dynamics(t, x, u(t)), time, x0);

%% Plot the Results
figure;
subplot(2,1,1);
plot(t, x(:, 1), 'r', t, x(:, 2), 'b');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('\theta_1', '\theta_2');
title('Joint Angles');

subplot(2,1,2);
plot(t, x(:, 3), 'r', t, x(:, 4), 'b');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('\dot{\theta}_1', '\dot{\theta}_2');
title('Joint Angular Velocities');
