function [C_w,H,H_zz,H_tt] = FRF(omega_1,U,B,rho,M,K,C)

C_w = C_omega(omega_1,U,B); % Fct circulatoire de Theodorsen approximée par Jones

q = pi*rho*U^2*B;
A = [B*omega_1^2/(4*U^2) - C_w*1i*omega_1/U , B/(4*U)*1i*omega_1 + C_w*(1+B/(4*U)*1i*omega_1) ; -B/4*C_w*1i*omega_1/U , B/4*(B^2*omega_1^2/(32*U^2) - B/(4*U)*1i*omega_1 +  C_w*(1+B/(4*U)*1i*omega_1))];     
H = inv(-omega_1^2*M + 1i*omega_1*C + K - q*A); % Fct de transfert
% H = X/fb c'est pour cela que fb n'est pas pris en compte dans notre fct de transfert;
H_zz = H(1,1);
H_tt = H(2,2);

end