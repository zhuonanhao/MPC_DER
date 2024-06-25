function x_new = rk4_symbolic(dynamics, h, x_old, u)
    k1 = dynamics(x_old, u);
    k2 = dynamics(x_old + h/2 * k1, u);
    k3 = dynamics(x_old + h/2 * k2, u);
    k4 = dynamics(x_old + h * k3, u);
    x_new = x_old + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

