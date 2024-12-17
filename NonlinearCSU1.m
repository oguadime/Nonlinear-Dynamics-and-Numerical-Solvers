% This is the NonlinearBEU code that shows monotonocity (U^{n+1}<=U^n) and
% energy stability (F(U^{n})<= F(U^{n+1}))
function NonlinearCSU1(u0, tau)
    if nargin < 2
        u0 = 2;
        tau = 0.9;
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

    % Function for the energy F(U), the integral of f(u)
    F = @(u) -(u.^2) / 2 + (u.^4) / 4 + 1/4;

    % Initialize energy array
    energy = zeros(1, Nsteps);
    energy(1) = F(u0); % Initial energy

    % Main body of code
    for i = 2:Nsteps
        f = @(unext) unext - unew(i-1) - tau*unew(i-1) + tau*(unext).^3;
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
    plot(t_ode45, u_ode45, 'o-', 'DisplayName', 'u_{ode45}');
    hold on;
    plot(t, unew, '*-', 'DisplayName', 'u_{CS}');
    xlabel('t');
    ylabel('u(t)');
    legend('show');
    title('Comparison of Convexity Splitting and ode45 Solutions');
    grid on;

    % Plot the energy function over time
    figure(2);
    plot(t, energy, 'r*-');
    xlabel('t');
    ylabel('F(U^n)');
    title('Energy Function F(U^n) over Time');
    grid on;


end
