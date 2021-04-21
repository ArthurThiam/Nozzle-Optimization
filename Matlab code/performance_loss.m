function loss_factor = performance_loss(c, d, dz)

a = -1.3360;
b =  2.3935;
%d = 0.0454;

R = @(z) (a + (sqrt(b+1000*c*z)/d))*0.001;

% Determine z_min (z_value for which the graphite-zirconium transition
% radius is reached)
z_min = 0;
r = 0;
R_transition = 72.15/1000.;

while abs(r-R_transition) > 0.01
    r = R(z_min);
    z_min = z_min + (0.01/1000.);
end

z_max = z_min + dz;

% Determine the nozzle exit angle
delta_z = 0.005;
slope = (R(z_max)-R(z_max-delta_z))/(delta_z);
exit_angle = atan(slope);

loss_factor = 1-(1-cos(exit_angle))/2;