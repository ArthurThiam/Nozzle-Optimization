function minimum = gradient_descent_simple(func, x0)
% Optimization parameters
max_iterations = 10000;
threshold = 1e-5;

% Optimization loop
i = 0;                  % Iteration counter
x = x0;                 % Initial design point
ffd = 1e-8;             % finite differences stepsize
l_rate = 0.0000001;     % learning rate
run = 1;

while (run == 1)
    
     % evaluate current function value
     f = func(x);
    
     % evaluate function at points left and right of previous point
     x_a = [x(1)+ffd x(2) x(3)];
     f_a = func(x_a);
     
     x_b = [x(1) x(2)+ffd x(3)];
     f_b = func(x_b);
     
     % evaluate gradients and descent direction
     grad_a = (f_a-func(x))/ffd;
     grad_b = (f_b-func(x))/ffd;
     
     grad = [grad_a grad_b];
     descent = -grad;
     fprintf('Steepest descent = ')
     disp(descent)
     
     step(1) = l_rate*descent(1);
     step(2) = l_rate*descent(2);
     fprintf('Step sizes = ')
     disp(step)
     
     % new design point
     x_new(1) = x(1) + step(1);
     x_new(2) = x(2) + step(2);
     x_new(3) = x(3);
     
     f_old = f;
     f = func(x_new);
     difference = f - f_old;
     fprintf('Objective value change = ')
     disp(difference);
     
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
     
     %plot(x(1), x(2), 'x', x_new(1), x_new(2), 'o')%,x_new(1),x_new(2),'o')
     x = [x_new x(3)];
     
     fprintf('previous function value = ')
     disp(f_old)
     fprintf('descent direction = ')
     disp(descent)
     fprintf('new function value = ')
     disp(f)
     fprintf('iteration = ')
     disp(i)
     fprintf('---------------------\n\n')
     
     %run = input('  Another optimization cycle? (0: No  1: Yes):');    
     %pause(0.5);
end   
     
minimum = x;