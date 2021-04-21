function U_e = exhaust_velocity(R_e, gamma, c, d, dz)

throat_area = pi*(40*0.001)^2;
R = 640.91;
T_c = 2170;

Gamma = sqrt(gamma)*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)));

expansion_ratio = (pi*R_e^2)/(throat_area);
p_ratio = pressure_ratio(expansion_ratio, gamma);

U_e = sqrt(2*(gamma*R*T_c)/(gamma-1)*(1-(p_ratio)^((gamma-1)/gamma)))*performance_loss(c, dz, d);