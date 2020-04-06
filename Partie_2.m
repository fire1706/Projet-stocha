%% Partie 2 => Analyse temporelle du problème
clear all;
%% DONNEES
%Propriétés du PONT:
k = 1;
m = 22740;
J = 2.47*10^6;
B = 31;
xhi = 0.003;
f_z = 0.10;
f_theta = 0.278;

%Propriétés du VENT:
rho = 1.22;
L_w = 20;
U = 10:1:76; % Vecteur utile => U = 10:1:76;  ATTENTION => Si vous générez une boucle de U, faites attention aux nombreuses figures qui seront plottées!!! 
I_w = 0.05;

%Valeurs utiles trouvées à la Partie 1:
w_0 = [2*pi*f_z,2*pi*f_z;2*pi*f_theta,2*pi*f_theta];
M = [m,0;0,J];
K = w_0.^2.*M;
C = 2*xhi*sqrt(K*M);% Egalement valable: C = 2*xhi*w_0.*M;

%Définition des variables pour cette partie:
t_max = 600;
N = t_max; % Nombre d'échantillons % Début : N = 2^(12) ???
time = linspace(0,t_max,N); % Vecteur temps

d_omega = 2*pi/t_max; 
omega_max = (N-1)*d_omega;
omega = (0:N-1)*d_omega;

%% Résolution

%A)Génération de réalisations du processus w:

j = 1;
while j <= length(U)
    
k = 1;
while k <= length(omega)
    theta(k) = 2*pi*rand(1,1);
    Mat_S(k) = S_omega(omega(k),U(j)); % PSD en fréquentiel 
    W(k) = N*sqrt(Mat_S(k)*d_omega)*exp(1i*theta(k)); % W(omega)
    k = k + 1;
end
W(1) = 0; % Conditions (pour vérifier la symétrie) pour la décomposition en série de Fourier d’une fonction continue: N est pair
W(N/2 + 1) = 0;
W(N/2 + 2 : N)  = conj(W(N/2:-1:2));
w = ifft(W,'symmetric');
% figure;
% plot(time,w,'LineWidth',1.2);
% xlabel('Temps [sec]');
% ylabel('Composante fluctuante du vent w');
% title('Processus aléatoire w(t)');

%B) Implémentation d’un schéma d’intégration pour résoudre l’équation du mouvement et calcul de la réponse par ODE45:
F_b_omega = 1/4*rho*U(j)*B*[4*pi;pi*B]*w;
L_b = F_b_omega(1,:);
M_b = F_b_omega(2,:);
sigma_w = std(w); % Ecart-type de w(t)
sigma_Sw = sqrt( trapz(omega,Mat_S)*2); % Ecart_type de la PSD en fréquentiel

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

p = 1;
while p <= length(time)
f_b(:,p) = [0;0;L_b(p);M_b(p);0;0;0;0;0;0]; % Forces de turbulence en temporelle 
p = p + 1;
end

% var = 2; % On démarre à var = 2 pour éviter de calculer la fct de transfert quand omega = 0 (sinon on a une matrice remplie de valeurs infinies) 
% while var <= length(omega)
% H(:,:,var-1) = inv(1i*omega(var)*M_final - A_final); % FCT de transfert
% H_zz(var-1) = H(1,1,var-1);
% H_tt(var-1) = H(2,2,var-1);
% var = var + 1;
% end
% valeurs_propres = real(eig(A_final)); % Pour la vérification que les parties réelles des valeurs propres soient négatives

% figure;
% plot(omega,abs(H_zz),'LineWidth',1.5); % H_zz
% xlim([0 2]);
% title( 'Nodal FRF H_{zz}(\omega)' )
% xlabel('Pulsations \omega [rad/s]')
% ylabel('FRF')
% figure;
% plot(omega,abs(H_tt),'LineWidth',1.5); % H_thetatheta
% xlim([0 2]);
% title( 'Nodal FRF H_{\theta\theta}(\omega)' )
% xlabel('Pulsations \omega [rad/s]')
% ylabel('FRF')

tspan = [0 t_max];
z0 = zeros(10,1); % structure initialement immobile
[t,z] = ode45(@(t,z) differentielle(t,z,A_final,M_final,f_b) ,tspan,z0); %ODE45 résout notre équation différentielle " M*z_dot = A*z + f_b"
% figure;
% plot(t,z(:,1),t,z(:,2),'LineWidth',1.5); % Graphe de la réponse z(t) et theta(t) de l'équation du mouvement
% xlabel('Temps [sec]');
% title("Réponse de l'équation du mouvement par ODE45");
% legend('z(t)','\theta(t)');

%C) Calcul de la variance de la réponse:
dt = time(2)-time(1);

[PSD_z,F_z] = pwelch(z(:,1),[],[],length(omega),1/dt); % PSD de z(t)
F_zz = 2*pi*F_z; % Transformations de temporel vers fréquentiel pour la comparaison dans le domaine fréquentiel
PSD_zz = PSD_z/(4*pi);

[PSD_t,F_t] = pwelch(z(:,2),[],[],length(omega),1/dt); % PSD de theta(t)
F_tt = 2*pi*F_t; % Transformations de temporel vers fréquentiel pour la comparaison dans le domaine fréquentiel
PSD_tt = PSD_t/(4*pi);

% figure;
% semilogy(F_zz,PSD_zz,'LineWidth',1.5);
% xlim([0 2]);
% xlabel('Pulsation \omega [rad/s]');
% ylabel('PSD');
% title(' S_{zz}(\omega)');
% figure;
% semilogy(F_tt,PSD_tt,'LineWidth',1.5);
% xlim([0 2]);
% xlabel('Pulsation \omega [rad/s]');
% ylabel('PSD');
% title(' S_{\theta\theta}(\omega)');

% figure;
% semilogy(F_zz,PSD_zz,omega,Mat_S,'LineWidth',1.5); % Comparaison entre l'approche temporelle et l'approche fréquentielle
% xlim([0 2]);
% xlabel('Pulsation \omega [rad/s]');
% ylabel('PSD');
% title(' S_w(\omega)');
% legend('Approche temporelle','Approche fréquentielle')

sigma_z(j) = std(z(:,1)); % Ecart-type de z(t)
sigma_theta(j) = std(z(:,2)); % Ecart-type de theta(t)

j = j + 1;
end

%D) Calcul de la vitesse de vent critique au flottement:

U = 10:1:100;
%U_crit = 76.86; % Valeur calculée dans la partie 1
sigma_z(68:length(U)) = 0; % Etant donné que le système devient instable passé la vitesse critique moyenne du vent U_crit, on va définir à zéro les valeurs d'écart-type pour U > U_crit
sigma_theta(68:length(U)) = 0;
figure;
plot(U,sigma_z,'-o',U,sigma_theta*B/2,'-o','LineWidth',1.5);
xlabel('Vitesse moyenne du vent U [m/s]');
ylabel('Ecart-type');
title('Ecart-type de la réponse');
legend('\sigma_{zz}','\sigma_{\theta\theta} * B/2')
