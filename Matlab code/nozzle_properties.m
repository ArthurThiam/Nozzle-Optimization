function [Result_2] = nozzle_properties(c, d, dz)

% Initialize parameters
%dz = 100/1000;% m

a = -1.3360;
b =  2.3935;
%c is called from the function input
%d =  0.0454;

% reverse build function
R = @(z) (a + (sqrt(b+1000*c*z)/d))*0.001;

% Determine z_min (z_value for which the graphite-zirconium transition
% radius is reached)
z_min = 0;
r = 0;
R_transition = 72.15/1000.;

while abs(r-R_transition) > 0.1/1000.
    r = R(z_min);
    z_min = z_min + 0.001/1000.;
end
    
z_max = z_min + dz;       % m

% Determine Results
V_zirconium = (pi*1./1000.*(1./1000.*dz + 2*integral(R, z_min, z_max)))*1000000; %cm3
V_titanium = (pi*4./1000.*(4./1000.*dz + 2*integral(R, z_min, z_max)))*1000000; %cm3
R_e = R(z_max);

% Return Results
Result_2 = [V_titanium V_zirconium R_e];