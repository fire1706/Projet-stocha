%% Partie 1 

%DONNEES:
%1) Propriétés PONT
i = 1;
w = 0:0.1:2;
m = 22740;
J = 2.47*10^6;
B = 31;
xhi = 0.003;
f_z = 0.10;
f_theta = 0.278;

%2) Propriétés VENT
rho = 1.22;
L_w = 20;
U = 10;%U = [0,100];
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
j = 1;
H = zeros(2,2);
while i <= length(w)    
    C_w(i) = C_omega(w(i),U,B);
    q = pi*rho*U^2*B;
    Grosse_matrice_desesmorts = [B*w(i)^2/(4*U^2) - C_w*1i*w(i)/U , B/(4*U)*1i*w(i) + C_w*(1+B/(4*U)*1i*w(i)) ; -B/4*C_w*1i*w(i)/U , B/4*(B^2*w(i)^2/(32*U^2) - B/(4*U)*1i*w(i) +  C_w*(1+B/(4*U)*1i*w(i)))];
    H(i,j) = (-w(i)^2*M + 1i*w(i)*C + K - q*Grosse_matrice_desesmorts).^(-1);
    j = j + 1;
% F_se_omega = q*Grosse_matrice_desesmorts*[POINT_interro1;POINT_interro2]; 
% F_b_omega = 1/4*rho*U*B*[4*pi;pi*B]*W_omega;
% H = X/fb;

%B) Calcul densité de puissance 
    S_w(i) = S_omega(w(i),U);
    i = i + 1;
end
%C)
M_conj = [pi*rho*B^2/4,0;0,pi*rho*B^4/128];
C_conj = [0,-q*B/(4*U);0,q*B^2/(16*U)];



