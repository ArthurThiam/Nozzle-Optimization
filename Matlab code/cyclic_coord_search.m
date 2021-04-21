function minimum = cyclic_coord_search(func, x0, x_l, x_u)
% Optimization parameters
max_iterations = 10000;
threshold = 1e-15;

% Optimization loop
i = 0;                  % Iteration counter
x = x0;                 % Initial design point
ffd = 1e-8;             % finite differences stepsize
l_rate = 0.0000001;     % learning rate
run = 1;

% Boundaries

while (run == 1)
    
     % evaluate current function value
     f = func(x);
     
    
     % set boundaries based on cycle.
     if (mod(i,3) == 0) 
        x_lower = [x_l(1) x(2) x(3)];
        x_upper = [x_u(1) x(2) x(3)];
     elseif (mod(i, 3) == 1)
        x_lower = [x(1) x_u(2) x(3)];
        x_upper = [x(1) x_l(2) x(3)];
     elseif (mod(i, 3) == 2)
        x_lower = [x(1) x(2) x_l(3)];
        x_upper = [x(1) x(2) x_u(3)];
     end
     
     % perform line search
     x_new = fminbnd(func, x_lower, x_upper);
     
     % Determine new value
     fprintf('New design point = ')
     disp(x_new);
     f_new = func(x_new);
     fprintf('New function value = ')
     disp(f_new);
     
     difference = f_new - f;
     fprintf('Objective value change = ')
     disp(difference);
     fprintf('Apogee = ')
     disp(1e4/f_new);   %km
     
     % verify termination tolerance
     if abs(difference) < threshold
         run = 0;
         fprintf('threshold value reached.\n')
     end
     
     % verify iteration count
     i = i + 1;
     if i >= max_iterations
         run = 0;
         fprintf('iteration count reached.\n')
     end
     
     % Update design point
     x = x_new;
     pause(1)
     fprintf('-----------------------------------\n');
end

fprintf('Optimum = ')
disp(x);
minimum = x;