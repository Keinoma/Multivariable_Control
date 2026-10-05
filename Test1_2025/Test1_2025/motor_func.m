function dx = motor_func(u, La, Ra, Km, Re, J, bm)

    % from the input vector u, we extract the states
    x(1) = u(1);
    x(2) = u(2);
    x(3) = u(3);

    % from the input vector u, we extract the control inputs of the system
    uin(1) = u(4);
    uin(2) = u(5);

    % we output the time derivative of the states
    dx(1) = 1/La*(uin(1) - Ra*x(1) - Km*x(2)*x(3));
    dx(2) = uin(2) - Re*7*atanh(x(2)/550);
    dx(3) = 1/J*(Km*x(1)*x(2) - bm*x(3));

end

