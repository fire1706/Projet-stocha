%% Partie 2 => Temporelle
%DONNEES:
%1) Propriétés PONT
k = 1;
%w = 0:0.01:2;
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
j = 1;
I_w = 0.05;

w_0 = [2*pi*f_z,2*pi*f_z;2*pi*f_theta,2*pi*f_theta];
M = [m,0;0,J];
K = w_0.^2.*M;
C = 2*xhi*sqrt(K*M);%C = 2*xhi*w_0.*M;

t_max = 600;
%Tps = 0:1:t_max;
N = t_max;
t = linspace(0,t_max,N);
d_omega = 2*pi/t_max;
omega_max = (N-1)*d_omega;
omega = 0:d_omega:omega_max;

k = 1;
%theta = 2*pi*rand(1,N);
while k <= length(omega)
    theta(k) = 2*pi*rand(1,1);
    S(k) = S_omega(omega(k),U);
    W(k) = N*sqrt(S(k)*d_omega)*exp(1i*theta(k));
    k = k + 1;
end
W(1) = 0;
W(N/2 + 1) = 0;
W(N/2 + 2 : N) = conj(W(N/2:-1:2));

w = ifft(W,'symmetric');
plot(t,w,'LineWidth',1);

F_b_omega = 1/4*rho*U*B*[4*pi;pi*B]*W;
L_b = F_b_omega(1,:);
M_b = F_b_omega(2,:);
ecart_type = std(W);
lol = std(S);
sigma_from_S_w = sqrt( trapz(omega,S) );
% %Vérif :
% a = sum(w.^2);
% b = 1/N * sum(abs(W).^2);
% if (a == b) 
%     True = 1;
% else 
%     True = 0;
% end


% M_conj = [pi*rho*B^2/4,0;0,pi*rho*B^4/128];
% C_conj = [0,-q*B/(4*U(j));0,q*B^2/(16*U(j))];
% I = [1,0;0,1];
% M_0 = [I,0;0,M+M_conj];
% A_0 = [0,I;-K,-(C+C_conj)];
% M = [M_0,0;-A,I];
% A = [A_0,S;0,B];
