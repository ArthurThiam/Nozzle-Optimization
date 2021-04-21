function atm = atmosphere(x)

if (0 <= x) && (x < 11000)
    x_b = 0;
    lapse = -0.0065;
    T_b = 288.15;
    rho_b = 1.2250;
    p_b = 101325;
    
elseif (11000 <= x) && (x < 20000)
    x_b = 11000;
    lapse = 0;
    T_b = 216.65;
    rho_b = 0.36391;
    p_b = 22632.1;
    
elseif (20000 <= x) && (x < 32000)
    x_b = 20000;
    lapse = 0.001;
    T_b = 216.65;
    rho_b = 0.08803;
    p_b = 5474.89;
    
elseif (32000 <= x) && (x < 47000)
    x_b = 32000;
    lapse = 0.0028;
    T_b = 228.65;
    rho_b = 0.01322;
    p_b = 868.02;
    
elseif (47000 <= x) && (x < 51000)
    x_b = 47000;
    lapse = 0;
    T_b = 270.65;
    rho_b = 0.00143;
    p_b = 110.91;
    
elseif (51000 <= x) && (x < 71000)
    x_b = 51000;
    lapse = -0.0028;
    T_b = 270.65;
    rho_b = 0.00086;
    p_b = 66.94;
    
else
    x_b = 71000;
    lapse = -0.002;
    T_b = 214.65;
    rho_b = 0.000064;
    p_b = 3.96;
    
end

if lapse == 0
    rho = rho_b*exp((-9.80665*0.0289644*(x-x_b))/(8.3144598*T_b));
    p = p_b*exp((-9.80665*0.0289644*(x-x_b))/(8.3144598*T_b));
    
elseif lapse ~= 0
    rho = rho_b*(T_b/(T_b+lapse*(x-x_b)))^(1+(9.80665*0.0289644)/(8.3144598*lapse));
    p = p_b*(T_b/(T_b+lapse*(x-x_b)))^((9.80665*0.0289644)/(8.3144598*lapse));
    
end

atm = [p rho];
