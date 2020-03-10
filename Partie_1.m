%% Partie 1 

%DONNEES:
%1) Propriétés PONT
m = 22740;
J = 2.47*10^6;
B = 31;
epsilon = 0.003;
f_z = 0.10;
f_theta = 0.278;

%2) Propriétés VENT
rho = 1.22;
L_w = 20;
%U = [0,100];
I_w = 0.05;

%A) Fct de transfert:
M = [m;J];

C = 2*epsilon*sqrt(K*M);

