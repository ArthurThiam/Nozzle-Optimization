function performance_cost = objective_function(x)         % Calculates the performance/cost ratio to be minimized
%% Initialize Parameters
vol_cost_titanium = 1;                                       % [€/m3] Titanium volumetric cost weight
vol_cost_zirconium = 3;                                      % [€/m3] Zirconium volumetric cost weight

%% Determine Cost
vol_list = volume(x(1), x(2), x(3));
cost = 10000*(vol_list(1)*vol_cost_titanium + vol_list(2)*vol_cost_zirconium);      % Total cost of the skirt
%fprintf('cost determined.\n');

%% Determine Performance
I_sp = specific_impulse(vol_list(3), 1.15, x(1), x(2), x(3));
skirt_mass = vol_list(1)*4506 + vol_list(2)*6520;
delta_V = 9.80665*I_sp*log((350+skirt_mass)/(350+skirt_mass-174-40));
%fprintf('performance determined.\n');

%% Print Results
%fprintf('c = ');
%disp(x(1));
%fprintf('d = ');
%disp(x(3));
%fprintf('dz = ');
%disp(x(2));
%fprintf('\n');

%fprintf('skirt_mass [kg] = ');
%disp(skirt_mass);
%fprintf('delta V [m/s] = ');
%disp(delta_V);
%fprintf('cost [€] = ');
%disp(cost);
%fprintf('----------------------------------\n');

%% Return Performance Cost
performance_cost(1) = cost;
performance_cost(2) = delta_V*10^5;
