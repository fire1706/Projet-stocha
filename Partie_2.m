%% Partie 2 => Temporelle
clear all;
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
N = t_max; % Début : N = 2^12 ???
time = linspace(0,t_max,N);
d_omega = 2/599;%2*pi/t_max;
omega_max = (N-1)*d_omega;
%omega = 0:d_omega:omega_max;
omega = (0:N-1)*d_omega;

while j <= length(U)
    
k = 1;
while k <= length(omega)
    theta(k) = 2*pi*rand(1,1);
    Mat_S(k) = S_omega(omega(k),U(j));
    W(k) = N*sqrt(Mat_S(k)*d_omega)*exp(1i*theta(k));
    k = k + 1;
end
figure;
semilogy(omega,Mat_S,'LineWidth',1.5);

w = ifft(W,'symmetric');
figure;
plot(time,w,'LineWidth',1.2);

F_b_omega = 1/4*rho*U(j)*B*[4*pi;pi*B]*w;
L_b = F_b_omega(1,:);
M_b = F_b_omega(2,:);
ecart_type = std(w); % ecart-type
sigma_from_S_w = sqrt( trapz(omega,Mat_S)*2); %integrale de la PSD en fréquentiel

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
%f_b = zeros;
while k <= length(omega)
f_b(:,k) = [0;0;L_b(k);M_b(k);0;0;0;0;0;0];
k = k + 1;
end
fb = interp1(time,f_b.',t).';
tspan = [0 t_max];
z0 = zeros(10,1); % structure initialement immobile
[t,z] = ode45(@(t,z) differentielle(t,z,A_final,M_final,f_b) ,tspan,z0);

% y(1,:) = transpose(z(:,1)); % composante z
% y(2,:) = transpose(z(:,2)); % composante theta
% y(3,:) = transpose(z(:,3));
% y(4,:) = transpose(z(:,4));
% v(1,:) = transpose(z(:,5));
% v(2,:) = transpose(z(:,6));
% v(3,:) = transpose(z(:,7));
% v(4,:) = transpose(z(:,8));
% v(5,:) = transpose(z(:,9));
% v(6,:) = transpose(z(:,10));

figure;
plot(t,z(:,1),t,z(:,2),'LineWidth',1.5);
ect_z(j) = std(z(:,1));
ect_theta(j) = std(z(:,2));

dt = time(2)-time(1);
[PSD_z,F_z] = pwelch(z(:,1),[],[],length(omega),1/dt);
[PSD_t,F_t] = pwelch(z(:,2),[],[],length(omega),1/dt);

F_zz = 2*pi*F_z;
PSD_zz = PSD_z/(4*pi);
F_tt = 2*pi*F_t;
PSD_tt = PSD_t/(4*pi);

figure;
hold on;
semilogy(F_zz,PSD_zz,'LineWidth',1.5);
semilogy(omega,Mat_S,'LineWidth',1.5);
hold off;

figure;
semilogy(F_tt,PSD_tt,'LineWidth',1.5);


j = j + 1;
end
figure;
plot(U,ect_z,'-o',U,ect_theta*B/2,'-o','LineWidth',1.5);

