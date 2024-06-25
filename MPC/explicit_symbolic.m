function x_new = explicit_symbolic(dynamics, h, x_old, u)
    x_new = x_old + h * dynamics(x_old, u);
end