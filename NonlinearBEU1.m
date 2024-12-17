% This is the NonlinearBEU code that shows monotonocity (U^{n+1}<=U^n) and
% energy stability (F(U^{n})<= F(U^{n+1})) and plots the energt function
function NonlinearBEU1(u0, tau)
    if nargin < 2
        u0 = 2;
        tau = 0.1;
    end    
    
    % Discretization Parameters
    T = 5;
    dt = tau;
    t = 0:dt:T;
    Nsteps = length(t);

    % Initial conditions
    uold = u0;
    unew = 0*t;
    unew(1) = uold;

    % Function for the energy F(U)
    F = @(u) -(u.^2) / 2 + (u.^4) / 4 + 1/4;

    % Initialize energy array
    energy = zeros(1, Nsteps);
    energy(1) = F(u0); % Initial energy

    % Main body of code
    for i = 2:Nsteps
        f = @(unext) unext - unew(i-1) - tau*(unext - unext^3);
        unew(i) = fzero(f, unew(i-1));
        
        % Calculate F(U^n)
        energy(i) = F(unew(i));
        
        % Display U^n and F(U^n)
        fprintf('Time step %d: U^n = %.4f, F(U^n) = %.4f\n', i-1, unew(i), energy(i));
    end   

    % Exact solution with Ode45
    [t_ode45, u_ode45] = ode45(@(t, u) u - u^3, [0 T], uold);

    % Plot the results
    figure(1);
    plot(t_ode45, u_ode45, 'o-', t, unew, '*-');
    xlabel('t');
    ylabel('u(t)');
    legend('u_{ode45}', 'u_{BE}');
    title('Comparison of Backward Euler and ode45 Solutions');
    grid on;

    % Plot the energy function over time
    figure(2);
    plot(t, energy, 'r*-');
    xlabel('t');
    ylabel('F(U^n)');
    title('Energy Function F(U^n) over Time');
    grid on;

end
