clear;clc;
import casadi.*

%% Define system
num_states = 2;
num_controls = 1;

x = SX.sym('x', num_states);
u = SX.sym('u', num_controls);

f = [(1-x(2)^2)*x(1) - x(2);
     x(1) + u];
df_dx = jacobian(f, x)