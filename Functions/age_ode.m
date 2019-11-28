function age_prime = age_ode(x, a, q, phi, Grid)
    % interpolates data sets at x
    q = interp1(Grid.xf, q, x, 'pchip', 'extrap');
    phi = interp1(Grid.xc, phi, x, 'pchip', 'extrap');
    age_prime = phi/q;
end