%% Partie 1 

%DONNEES:
%1) Propriétés PONT
m = 22740;
J = 2.47*10^6;
B = 31;
xhi = 0.003;
f_z = 0.10;
f_theta = 0.278;

%2) Propriétés VENT
rho = 1.22;
L_w = 20;
%U = [0,100];
I_w = 0.05;

%A) Fct de transfert:
w_0 = [2*pi*f_z; 2*pi*f_theta];
M = [m,0;0,J];
K = w_0.^2.*M; 
C = 2*xhi*sqrt(K*M);%C = 2*xhi*w_0.*M;

% M*x_ddot + C*x_dot + K*x = f
a0 = 1;
a1 = -0.165;
a2 = -0.335;
b1 = 0.0455;
b2 = 0.3;
w_1 = 2*U*b1/B;
w_2 = 2*U*b2/B;
C_omega = a0 + a1*(1i*omega)/(1i*omega + w_1) + a2*(1i*omega)/(1i*omega + w_2);
q = pi*rho*U^2*B;
F_se_omega = q*[B*omega^2/(4*U^2) - C_omega*1i*omega/U , B/(4*U)*1i*omega + C_omega*(1+B/(4*U)*1i*omega) ; -B/4*C_omega*1i*omega/U , B/4*(B^2*omega^2/(32*U^2) - B/(4*U)*1i*omega +  C_omega*(1+B/(4*U)*1i*omega))]*[POINT_interro1;POINT_interro2]; 
F_b_omega = 1/4*rho*U*B*[4*pi;pi*B]*W_omega;

