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
% A = [1 1;0 1];
% B = [0; 1];
% rhs = A * x + B * u;

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

X(:,1) = X0;
cost = 0;
for k = 1:N
    X(:,k+1) = X(:,k) + dt * F(X(:,k), U(:,k));
    cost = cost + (X(:,k)-Xf)'*Q*(X(:,k)-Xf) + U(:,k)'*R*U(:,k);
end

OPT_variables = reshape(U, num_controls * N, 1);

nlp = struct;
nlp.f = cost;
nlp.x = OPT_variables;
nlp.p = [X0; Xf];

opts = struct;
opts.ipopt.print_level = 0;
opts.print_time = 0;

solver = nlpsol('solver', 'ipopt', nlp, opts);

% Simulation params
dt = 0.1;

tf = 20;
time = 0:dt:tf;
x0 = [3;0];
xf = [1;0];
xold = x0;
figure;
hold on
for i = 2:length(time)
    p = [xold;xf];
    u0 = zeros(N,num_controls);
    sol = solver('x0',reshape(u0,num_controls*N,1),'p',p);
    u_array = reshape(full(sol.x)', num_controls, N);
    u = u_array(1);
    xnew = xold + dt * F(xold, u);
    xold = xnew;
    scatter(full(xold(1)),full(xold(2)))

end
