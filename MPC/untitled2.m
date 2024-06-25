clear;clc;

syms t u(t) x1(t) x2(t) lambda1(t) lambda2(t)
assume([t u(t), x1(t), x2(t), lambda1(t), lambda2(t)], 'real')

x = [x1 x2]';
lambda = [lambda1 lambda2]';

f = [(1-x2^2)*x1-x2;
     x1 + u];

L = x' * eye(2) * x + u^2;

dL_dx = jacobian(L, x);
dL_du = jacobian(L, u);
df_dx = jacobian(f, x);
df_du = jacobian(f, u);

dlambda_dt = -dL_dx + lambda'*df_dx
dJ_du = dL_du - lambda' * df_du

