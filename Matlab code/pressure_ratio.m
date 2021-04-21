function p_ratio = pressure_ratio(epsilon_true, gamma)

Gamma = sqrt(gamma)*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
pressure_ratio_calculated = 0.001;
pressure_stepsize = 0.00005;
epsilon_calculated = 0;

error = abs(epsilon_true - epsilon_calculated);
threshold = 1;

while error > threshold
    epsilon_calculated = Gamma/sqrt((2*gamma/(gamma-1))*pressure_ratio_calculated^(2/gamma)*(1-pressure_ratio_calculated^((gamma-1)/gamma)));
    error = abs(epsilon_true - epsilon_calculated);
    
    pressure_ratio_calculated = pressure_ratio_calculated + pressure_stepsize;
    
end

p_ratio = pressure_ratio_calculated;