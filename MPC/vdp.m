% Define symbolic variables
syms t x1 x2 u lambda1 lambda2

% Define the state vector and adjoint vector
x = [x1; x2];
lambda = [lambda1; lambda2];

% Define the function f
f = [(1 - x2^2)*x1 - x2;
     x1 + u];

% Define the cost function L
L = x1^2 + x2^2 + u^2;

% Compute the Jacobian of f with respect to x
df_dx = jacobian(f, x);

% Compute the Jacobian of f with respect to u
df_du = jacobian(f, u);

% Compute the partial derivatives of L with respect to x and u
dL_dx = jacobian(L, x).';
dL_du = jacobian(L, u);

% Define the adjoint equations: d(lambda)/dt = -dL/dx + (df/dx)^T * lambda
adjoint_eq = diff(lambda, t) == -dL_dx + df_dx.' * lambda;

% Compute the gradient of the cost function with respect to the control: dJ/du
dJ_du = dL_du - df_du.' * lambda;

% Convert symbolic expressions to MATLAB functions for numerical computation
f_fcn = matlabFunction(f, 'Vars', {t, x1, x2, u});
adjoint_eq_fcn = matlabFunction(adjoint_eq, 'Vars', {t, x1, x2, lambda1, lambda2, u});
dJ_du_fcn = matlabFunction(dJ_du, 'Vars', {t, x1, x2, lambda1, lambda2, u});

% Simulation parameters
T = 10; % End time
dt = 0.01; % Time step
N = T/dt; % Number of time steps

% Initial conditions
x1_0 = rand; % Random initial condition for x1
x2_0 = rand; % Random initial condition for x2

% Initialize state and control variables
x1 = zeros(1, N); x2 = zeros(1, N);
x1(1) = x1_0; x2(1) = x2_0;
u = zeros(1, N); % Initial guess for control

% Optimization parameters
max_iter = 100;
alpha = 0.01; % Learning rate

% Time vector
time = linspace(0, T, N);

for iter = 1:max_iter
    % Forward solve for state equations
    for k = 1:N-1
        dx = f_fcn(time(k), x1(k), x2(k), u(k));
        x1(k+1) = x1(k) + dt * dx(1);
        x2(k+1) = x2(k) + dt * dx(2);
    end
    
    % Backward solve for adjoint equations
    lambda1 = zeros(1, N); lambda2 = zeros(1, N);
    for k = N:-1:2
        d_lambda = adjoint_eq_fcn(time(k), x1(k), x2(k), lambda1(k), lambda2(k), u(k));
        lambda1(k-1) = lambda1(k) - dt * d_lambda(1);
        lambda2(k-1) = lambda2(k) - dt * d_lambda(2);
    end
    
    % Compute the gradient of the cost function with respect to the control
    grad_u = zeros(1, N);
    for k = 1:N
        grad_u(k) = dJ_du_fcn(time(k), x1(k), x2(k), lambda1(k), lambda2(k), u(k));
    end
    
    % Update control using gradient descent
    u = u - alpha * grad_u;
end

% Plot results
figure;
subplot(3,1,1);
plot(time, x1);
title('State x1(t)');
xlabel('Time t');
ylabel('x1(t)');

subplot(3,1,2);
plot(time, x2);
title('State x2(t)');
xlabel('Time t');
ylabel('x2(t)');

subplot(3,1,3);
plot(time, u);
title('Control u(t)');
xlabel('Time t');
ylabel('u(t)');
