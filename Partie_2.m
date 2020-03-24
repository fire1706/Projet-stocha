%% Partie 2 => Temporelle
%DONNEES:
%1) Propriétés PONT
i = 1;
w = 0:0.01:2;
m = 22740;
J = 2.47*10^6;
B = 31;
xhi = 0.003;
f_z = 0.10;
f_theta = 0.278;

%2) Propriétés VENT
rho = 1.22;
L_w = 20;
U = 0:10:100;%U = [0,100];
j = 1;
I_w = 0.05;

w_0 = [2*pi*f_z,2*pi*f_z;2*pi*f_theta,2*pi*f_theta];
M = [m,0;0,J];
K = w_0.^2.*M;
C = 2*xhi*sqrt(K*M);%C = 2*xhi*w_0.*M;


M_conj = [pi*rho*B^2/4,0;0,pi*rho*B^4/128];
C_conj = [0,-q*B/(4*U(j));0,q*B^2/(16*U(j))];
I = [1,0;0,1];
M_0 = [I,0;0,M+M_conj];
A_0 = [0,I;-K,-(C+C_conj)];
M = [M_0,0;-A,I];
A = [A_0,S;0,B];
