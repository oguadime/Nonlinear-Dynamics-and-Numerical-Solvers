% Define the function F(u) and its derivative F'(u)
F = @(u) u.^3 - u;
Fprime = @(u) 3*u.^2 - 1;

% Convergence parameters
tol = 1e-12; % Tolerance for convergence
maxit = 10;  % Maximum number of iterations

% Narrow range of initial guesses near u = 1
initial_guesses = 0.5:0.1:1.2;  % Refined range of initial guesses
converged_to_one = [];           % Array to store initial guesses that converge to u = 1

% Loop over each initial guess
for u0 = initial_guesses
    u = u0;  % Set the initial guess
    iter = 0;
    
    % Newton-Raphson iteration loop
    while 1
        % Evaluate the function and its derivative at the current guess
        Fu = F(u);
        Fpu = Fprime(u);
        
        % Update the guess using Newton's method
        u_new = u - Fu / Fpu;
        
        % Check for convergence
        if abs(u_new - u) < tol
            % Check if the final value is close to 1
            if abs(u_new - 1) < tol
                converged_to_one = [converged_to_one, u0]; % Store the initial guess
            end
            break; % Convergence achieved, exit the loop
        end
        
        % Update for next iteration
        u = u_new;
        iter = iter + 1;
        
        % Check for maximum iterations
        if iter > maxit
             error ('no convergence')
        end
    end
end

% Display the initial guesses that converge to u = 1
disp('Initial guesses that converge to u = 1:');
disp(converged_to_one);
