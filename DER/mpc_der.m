clear;clc;close all

import casadi.*
% 
% 
% % Symbols/expressions
% x = MX.sym('x');
% y = MX.sym('y');
% z = MX.sym('z');
% f = x^2+100*z^2;
% g = z+(1-x)^2-y;
% 
% nlp = struct;            % NLP declaration
% nlp.x = [x;y;z];         % decision vars
% nlp.f = f;               % objective
% nlp.g = g;               % constraints
% 
% opts = struct;
% opts.ipopt.print_level = 0; % Suppress IPOPT output
% opts.print_time = 0;        % Suppress CasADi timing info
% 
% F = nlpsol('F','ipopt',nlp,opts);
% 
% % Solve the problem using a guess
% r = F('x0',[2.5 3.0 0.75],'ubg',0,'lbg',0);
% x_opt = r.x
% 
% import casadi.*


x = MX.sym('x',2); % Two states

rhs = [(1-x(2)^2)*x(1)-x(2);
       x(1)];

ode = struct;    % ODE declaration
ode.x   = x;     % states
ode.ode = rhs;   % right-hand side

% Construct a Function that integrates over 4s
F = integrator('F','cvodes',ode,0,4);

% Start from x=[0;1]
res = F('x0',[0;1]);

disp(res.xf)

% Sensitivity wrt initial state
res = F('x0',x);
S = Function('S',{x},{jacobian(res.xf,x)});

disp(S([0;1]))