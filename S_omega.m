function sigma_w = S_omega(w)
U = 0:1:100;
L_w = 20;
I_w = 0.05;
i = 0;
S = zeros;
while i <= 100
    sigma_w = U(i)*I_w;
    S(i) = sigma_w^2*L_w/(pi*U(i)) *(1+755.2*(w*L_w/(2*pi*U(i)))^2)/(1+283.2*(w*L_w/(2*pi*U(i)))^2)^(11/6);
end
end


