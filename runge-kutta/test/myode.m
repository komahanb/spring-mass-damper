function dz = myode(v,z)
dz = zeros(2,1);    
dz(1) = exp(z(2)) + sin(z(3));   
dz(2) = -0.5*z(1) + 0.5*exp(z(3));
dz(3) = -0.5*exp(z(2)) + 2.5*z(3);
