%% Clear cache
clear; clc; close all

import casadi.*

%% Define system
num_states = 3;
num_controls = 1;

x = SX.sym('x', num_states);
u = SX.sym('u', num_controls);

% rhs = [(1-x(2)^2)*x(1) - x(2);
%         x(1) + u];

A = [1 1 0;
     0 -1 1;
     1 0 -1];
B = [0; 1; 0];
rhs = A * x + B * u;

F = Function('F', {x, u}, {rhs});

%% Define MPC
dt = 0.1;
N = 10;

Q = eye(num_states);
R = 0.1*eye(num_controls);
P = 2*Q;
solver = initilizeMPC(F,N,dt,Q,R,P);

%% Simulation 

tf = 10;
time = 0:dt:tf;
timesteps = length(time);

x0 = rand(num_states,1);
xold = x0;

figure;

subplot(1,2,1)
hold on
xlim([0,tf])

subplot(1,2,2)
hold on

for i = 1:timesteps
    ctime = time(i);

    subplot(1,2,1)
    scatter(ctime, full(xold(1)), 'r.')

    subplot(1,2,2)
    scatter(full(xold(1)), full(xold(2)), 'r.')
    
    xf = [cos(ctime);0;0];
    p = [xold;xf];
    u0 = zeros(N,num_controls);
    sol = solver('x0',reshape(u0,num_controls*N,1),'p',p);
    u_array = reshape(full(sol.x)', num_controls, N);
    u = u_array(1);
    xnew = full(rk4_symbolic(F,dt,xold,u));
    xold = xnew;
end
