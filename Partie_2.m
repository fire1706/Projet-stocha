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

while j <= length(U)
k = 1;
%theta = 2*pi*rand(1,N);
while k <= length(omega)
    theta(k) = 2*pi*rand(1,1);
    Mat_S(k) = S_omega(omega(k),U(j));
    W(k) = N*sqrt(Mat_S(k)*d_omega)*exp(1i*theta(k));
    k = k + 1;
end
W(1) = 0;
W(N/2 + 1) = 0;
W(N/2 + 2 : N) = conj(W(N/2:-1:2));

w = ifft(W,'symmetric');
plot(t,w,'LineWidth',1);

F_b_omega = 1/4*rho*U(j)*B*[4*pi;pi*B]*w;
L_b = F_b_omega(1,:);
M_b = F_b_omega(2,:);
ecart_type = std(w);
sigma_from_S_w = sqrt( trapz(omega,Mat_S)*2);

b0 = 0;
b1 = 0.0455;
b2 = 0.3;
Mat_B = -2*U(j)/B *diag([b0 b1 b2 b0 b1 b2]);

a0 = 1;
a1 = -0.165;
a2 = -0.335;

Mat_A =[0,a0,0,0;0,a1,0,0;0,a2,0,0;0,0,-a0/U(j),a0*B/(4*U(j));0,0,-a1/U(j),a1*B/(4*U(j));0,0,-a2/U(j),a2*B/(4*U(j))];

q = pi*rho*U(j)^2*B;
S_matrice = q*[0,0,0,0,0,0;0,0,0,0,0,0;1,1,1,1,1,1;B/4,B/4,B/4,B/4,B/4,B/4];

M_conj = [pi*rho*B^2/4,0;0,pi*rho*B^4/128];
C_conj = [0,-q*B/(4*U(j));0,q*B^2/(16*U(j))];

M_0 = [eye(2),zeros(2,2);zeros(2,2),M+M_conj];
A_0 = [zeros(2,2),eye(2);-K,-(C+C_conj)];

M_final = [M_0,zeros(4,6);-Mat_A,eye(6)];
A_final = [A_0,S_matrice;zeros(6,4),Mat_B];

k = 1;
while k <= length(omega)
f_b(:,k) = [0;0;L_b(k);M_b(k);0;0;0;0;0;0];
%f_b2(:,:,k) = [0;0;L_b(k);M_b(k);0;0;0;0;0;0];
k = k + 1;
end

tspan = [0 t_max];
z0 = zeros(10,1); % structure initialement immobile
[t,z] = ode45(@(t,z) differentielle(t,z,A_final,M_final,f_b) ,tspan,z0);

j = j + 1;
end

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
