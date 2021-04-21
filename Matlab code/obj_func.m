function performance = obj_func(x)
gamma = 1.15;
m = 10;  %kg/s
altitude = 0; %m
t = 0; %s
dt = 0.1;
t_simulation = 400;
t_thrust = 22.5;
p_c = 15*10^5;  %Pa
A_t = pi*(40./1000.)^2;
c_d = 1.5;
v = 0;  % m/s
mass = 350; %kg

%% Determine nozzle parameters from input
nozzle = nozzle_properties(x(1), x(2), x(3));
A_e = (pi*nozzle(3)^2);

u_e = exhaust_velocity(nozzle(3), gamma, x(1), x(2), x(3));
loss_factor = performance_loss(x(1), x(2), x(3));

epsilon = A_e/A_t;
p_e = p_c*pressure_ratio(epsilon, gamma);

data = [];
time = [];
%% Start simulation loop
while (t < t_simulation)
   % Determine current atmospheric conditions
   atm = atmosphere(altitude);
   p_a = atm(1);
   rho = atm(2);
   
   
   % Calculate net force and acceleration
   if (t < t_thrust)
       F_t = loss_factor*m*u_e + (p_e-p_a)*A_e;
   else
       F_t = 0;
   end
   F_d = c_d*0.5*rho*v^2*((100./1000.)^2*pi);
   F_w = mass*9.80665;
   
   F = F_t - F_d - F_w;
   a = F/mass;
   
   % Calculate distance travelled
   d = v*dt + 0.5*a*dt^2;
   altitude = altitude + d;
   
   if (altitude < 0)
       altitude = 0;
   end
   
   % Update variables
   if (t < t_thrust)
       mass = mass - m*dt;
   end
   
   t = t + dt;
   v = v + a*dt;
   time = [time t];
   data = [data altitude];
end

%plot(time,data)
%fprintf('performance evaluated.\n');
performance = (1/max(data))*1e6;