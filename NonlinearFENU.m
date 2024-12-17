function NonlinearFENU(u0, tau)
    if nargin < 2
        u0 = 2;
        tau = 0.9;
    end    
    
    % Discretization Parameters
    T = 5;              % End time
    tol = 1e-4;         % Tolerance for error control
    t = 0;              % Start time
    unew = u0;          % Initial condition
    tsteps = t;         % Array to store time steps
    uFE = unew;        % Array to store solution

    % Main loop with variable time stepping
    while t < T
        % Forward Euler step
        u_next = unew + tau * (unew - unew^3);
        
        % Estimate error using two half steps
        u_half = unew + 0.5 * tau * (unew - unew^3);
        u_half_step = u_half + 0.5 * tau * (u_half - u_half^3);
        
        % Calculate the error estimate
        error_estimate = abs(u_next - u_half_step);
        
        % Adjust the time step based on error
        if error_estimate > tol
            tau = tau * 0.5;  % Reduce time step size if the error is too large
        else
            % Accept the step
            t = t + tau;
            unew = u_next;
            tsteps(end+1) = t;  % Store the time step
            uFE(end+1) = unew; % Store the solution
            
            % Increase time step size if the error is significantly smaller
            if error_estimate < tol * 0.1
                tau = tau * 1.5;
            end
        end
        
        % Adjust the last step to exactly reach T
        if t + tau > T
            tau = T - t;
        end
    end    

    % Exact solution with ode45
    [t_ode45, u_ode45] = ode45(@(t, u) u - u^3, [0 T], u0);

    % Plot the results
    figure;
    plot(t_ode45, u_ode45, 'o-', 'DisplayName', 'u_{ode45}');
    hold on;
    plot(tsteps, uFE, '*-', 'DisplayName', 'u_{FE}');
    xlabel('t');
    ylabel('u(t)');
    legend('show');
    title('Comparison of Forward Euler (Variable Step) and ode45 Solutions');
    grid on;
end
