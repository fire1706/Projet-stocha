function [S_zz,S_tt] = PSD_omega_1(omega_1,U,B,rho,H,L_w,I_w)

S_w = S_omega(omega_1,U,L_w,I_w); % Calcul de la PSD de Von Karman
F_b_omega = (1/4)*rho*U*B*[4*pi;pi*B];
S_b = F_b_omega*S_w*conj(transpose(F_b_omega));
S_x = H*S_b*conj(transpose(H));
S_zz = S_x(1,1); % PSD_zz
S_tt = S_x(2,2); % PSD_\theta\theta

end