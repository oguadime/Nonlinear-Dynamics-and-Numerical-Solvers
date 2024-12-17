function NonlinearCSNU(u0, tau)
    if nargin < 2
        u0 = 2;
        tau = 0.1;
    end    
    
    % Discretization Parameters
    T = 5;              % End time
    tol = 1e-4;         % Tolerance for error control
    t = 0;              % Start time
    unew = u0;          % Initial condition
    tsteps = t;         % Array to store time steps
    uCS = unew;        % Array to store solution

    % Main loop with variable time stepping
    while t < T
        % Define the implicit function for Convexity Splitting
        f = @(unext) unext - unew - tau * (unew - unext^3);
        
        % Solve for the next time step using fzero
        u_next = fzero(f, unew);  % Use the previous value as the initial guess
        
        % Estimate error by performing two half steps
        f_half = @(u_half) u_half - unew - 0.5 * tau * (unew - u_half^3);
        u_half = fzero(f_half, unew);
        f_half_next = @(u_half_next) u_half_next - u_half - 0.5 * tau * (u_half - u_half_next^3);
        u_half_step = fzero(f_half_next, u_half);
        
        % Calculate error estimate
        error_estimate = abs(u_next - u_half_step);
        
        % Adjust the time step based on the error
        if error_estimate > tol
            tau = tau * 0.5;  % Reduce the time step size if the error is too large
        else
            % Accept the step
            t = t + tau;
            unew = u_next;
            tsteps(end+1) = t;  % Store the time step
            uCS(end+1) = unew; % Store the solution
            
            % Increase the time step size if the error is sufficiently small
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
    plot(tsteps, uCS, '*-', 'DisplayName', 'u_{CS}');
    xlabel('t');
    ylabel('u(t)');
    legend('show');
    title('Comparison of Convexity Splitting (Variable Step) and ode45 Solutions');
    grid on;
end
