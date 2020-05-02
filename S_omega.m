function S = S_omega(w,U,L_w,I_w)
%Échelle de turbulence [m]
L_w = 20;
%Intensité de turbulence [-]
I_w = 0.05;

sigma_w = U*I_w; % Ecart-type

S = sigma_w^2*L_w/(pi*U) *(1+755.2*(w*L_w/(2*pi*U))^2)/(1+283.2*(w*L_w/(2*pi*U))^2)^(11/6); % PSD de Von Karman

end


