clear;clc;close all

% Parameters
T = 10;           % Time horizon
dt = 0.1;        % Time step
N = T / dt;       % Number of time steps
x0 = [1 0];           % Initial condition for x
u = zeros(N, 1);  % Initial guess for control

% Optimization parameters
max_iter = 100;   % Maximum number of iterations
alpha = 0.01;     % Learning rate

for iter = 1:max_iter
    % Forward solve for state equation
    x = zeros(N, 2);
    x(1,:) = x0;
    for k = 1:N-1
        x(k+1,:) = x(k,:) + dt * [(1-x(k,2)^2)*x(k,1) - x(k,2), x(k,1) + u(k)];
    end

    % Backward solve for adjoint equation
    lambda = zeros(N, 2);
    for k = N:-1:2
        lambda(k-1,:) = lambda(k,:) - dt * [lambda(k,2) - 2*x(k,1) - lambda(k,1)*(x(k,2)^2 - 1), - 2*x(k,2) - lambda(k,1)*(2*x(k,1)*x(k,2) + 1)];
    end

    % Compute gradient of cost function
    dJdu = zeros(N, 1);
    for k = 1:N
        dJdu(k) = 2*u(k) - lambda(k,2);
    end

    % Update control using gradient descent
    u = u - alpha * dJdu;
end

% Plot results
t = 0:dt:T-dt;
figure;
subplot(2,1,1);
plot(t, x);
title('State x(t)');
xlabel('Time t');
ylabel('x(t)');

subplot(2,1,2);
plot(t, u);
title('Control u(t)');
xlabel('Time t');
ylabel('u(t)');

