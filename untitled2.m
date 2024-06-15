%% Clear cache
clear; clc; close all

%% Initialization
import casadi.*

dt = 0.1;

num_states = 2;
num_controls = 1;

x = SX.sym('x', num_states);
u = SX.sym('u', num_controls);

rhs = [(1-x(2)^2)*x(1) - x(2);
        x(1) + u];

% Dynamics definition: x_prime = F(x,u) = rhs
F = Function('F', {x, u}, {rhs});

N = 10; % Horizon length

X = SX.sym('X', num_states, N+1);
U = SX.sym('U', num_controls, N);

% Initial and final conditions
X0 = SX.sym('X0', num_states);
Xf = SX.sym('Xf', num_states);

Q = eye(num_states);
R = eye(num_controls);

% Define the cost and constraints
cost = 0;
constraints = [];

% Initial state constraint
constraints = [constraints; X(:,1) - X0];

for k = 1:N
    % Dynamics constraint
    constraints = [constraints; X(:,k+1) - (X(:,k) + dt * F(X(:,k), U(:,k)))];
    % Cost function
    cost = cost + (X(:,k)-Xf)'*Q*(X(:,k)-Xf) + U(:,k)'*R*U(:,k);
end
% Terminal cost
cost = cost + (X(:,N+1)-Xf)'*Q*(X(:,N+1)-Xf);

% Reshape variables for optimization
OPT_variables = [reshape(X, num_states*(N+1), 1); reshape(U, num_controls*N, 1)];

% Define the NLP problem
nlp = struct('f', cost, 'x', OPT_variables, 'g', constraints, 'p', [X0; Xf]);

opts = struct;
opts.ipopt.print_level = 0;
opts.print_time = 0;

solver = nlpsol('solver', 'ipopt', nlp, opts);

% Simulation params
dt = 0.1;

tf = 20;
time = 0:dt:tf;
x0 = [0;0];
xf = [1;0];
xold = x0;

figure;
hold on;
for i = 2:length(time)
    p = [xold; xf];
    u0 = zeros(N, num_controls);
    x0_guess = repmat(xold, 1, N+1);
    initial_guess = [x0_guess(:); u0(:)];
    
    sol = solver('x0', initial_guess, 'lbx', -inf, 'ubx', inf, 'lbg', 0, 'ubg', 0, 'p', p);
    opt = full(sol.x);
    
    x_opt = reshape(opt(1:num_states*(N+1)), num_states, N+1);
    u_opt = reshape(opt(num_states*(N+1)+1:end), num_controls, N);
    
    u = u_opt(1);
    xnew = xold + dt * full(F(xold, u));
    xold = xnew;
    
    scatter(full(xold(1)), full(xold(2)));
end

xlabel('x1');
ylabel('x2');
title('State Trajectory');
hold off;
