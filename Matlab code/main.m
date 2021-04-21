clear all;
close all;
clc;

%% Initial Values

c = 0.0888;%0.0909;                                % [-] c value for nozzle curve definition
d = 0.0694;%0.0454;
dz = 0.1015;%0.1;                                  % [m]   skirt length
x0 = [c d dz];


%% Optimization Parameters

x_lower = [0.02 0.02 90./1000.];
x_upper = [0.2 0.1 120./1000.];


%% Optimization

f_objective = @(x) obj_func(x);
%x = fmincon(f_objective,x0,[],[],[],[],x_lower,x_upper)
%x = fminsearch(f_objective, x0)

%opt = cyclic_coord_search(f_objective, x0, x_lower, x_upper);

%% Plotting

c = 0.03:0.001:0.06;
d = 0.03:0.001:0.06;

Z = obj_func([c d dz]);
size([X,Y])
Z
mesh(X, Y, Z);
grid on;
hold on;

i = 0.04;
ii = 1;
j = 0.03;
ji = 1;
row = 1;

while i <= 0.05
    while j <= 0.06        
        f = 1e6/f_objective([i j 0.1]);

    end
    i = i + 0.001;
    ii = ii + 1;
    j = 0.03;    
end

Z
length(Z);

[X, Y] = meshgrid(c, d);
length([X,Y])

surf(X, Y, Z);










    