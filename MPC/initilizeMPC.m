function MPC_solver = initilizeMPC(dynamics,N,dt,Q,R,P)

import casadi.*

num_states = size(Q, 1);
num_controls = size(R, 1);

X = SX.sym('X', num_states, N+1);
U = SX.sym('U', num_controls, N);

X0 = SX.sym('X0', num_states);
Xf = SX.sym('Xf', num_states);

X(:,1) = X0;
cost = 0;

for k = 1:N
    X(:,k+1) = rk4_symbolic(dynamics, dt, X(:,k), U(:,k));
    cost = cost + (X(:,k)-Xf)'*Q*(X(:,k)-Xf) + U(:,k)'*R*U(:,k);
end
cost  = cost + (X(:,N+1)-Xf)'*P*(X(:,N+1)-Xf);

OPT_variables = reshape(U, num_controls * N, 1);

nlp = struct;
nlp.f = cost;
nlp.x = OPT_variables;
nlp.p = [X0; Xf];

opts = struct;
opts.ipopt.print_level = 0;
opts.print_time = 0;

MPC_solver = nlpsol('solver', 'ipopt', nlp, opts);

end