%% Partie 2 => Analyse temporelle du probl�me
clear all;
%% DONNEES
%Propri�t�s structurelles du PONT:
%Masse par unit� de longueur [kg/m]
m = 22740;
%Moment d�inertie par unit� de longueur [kg.m^2/m]
J = 2.47*10^6;
%Largeur du tablier [m]
B = 31;
%Amortissement structurel [-]
xhi = 0.003;
%Fr�quence propre verticale [Hz]
f_z = 0.10;
%Fr�quence propre de torsion [Hz]
f_theta = 0.278;

%Propri�t�s du VENT incident:
%Masse volumique de l�air [kg/m^3]
rho = 1.22;
%�chelle de turbulence [m]
L_w = 20;
%Vitesse moyenne horizontale du vent [m/s]
U = 10:1:50; % Vecteur utile => U = 10:1:76;  ATTENTION => Si vous g�n�rez une boucle de U, faites attention aux nombreuses figures qui seront plott�es!!! 
%Intensit� de turbulence [-]
I_w = 0.05;

%Valeurs utiles trouv�es � la Partie 1:
%Pulsation propre du tablier de pont [rad/s]
w_0 = [2*pi*f_z,2*pi*f_z;2*pi*f_theta,2*pi*f_theta];
%Matrice de masse
M = [m,0;0,J];
%Matrice de raideur
K = w_0.^2.*M;
%Matrice d'amortissement
C = 2*xhi*sqrt(K*M);% Egalement valable: C = 2*xhi*w_0.*M;

%D�finition des variables pour cette partie:
%Temps maximal d�sir� [sec]
t_max = 600;
%Nombre d'�chantillons
N = 600; % D�but: N = 2^(12) 
%Vecteur temps [sec]
time = linspace(0,t_max,N);

%Pas fr�quentiel [rad/s]
d_omega = 2*pi/t_max;
%Pulsation maximale d�sir�e [rad/s]
omega_max = (N-1)*d_omega;
%Pulsation [rad/s]
omega = (0:N-1)*d_omega;
omega(1) = 1e-03;
%% R�solution

%A)G�n�ration de r�alisations du processus w:

j = 1;
while j <= length(U)
    
k = 1;
while k <= length(omega)
    theta(k) = 2*pi*rand(1,1);
    Mat_S(k) = S_omega(omega(k),U(j)); % PSD en fr�quentiel 
    W(k) = N*sqrt(Mat_S(k)*d_omega)*exp(1i*theta(k)); % W(omega)
    k = k + 1;
end
W(1) = 0; % Conditions (pour v�rifier la sym�trie) pour la d�composition en s�rie de Fourier d�une fonction continue: N est pair
W(N/2 + 1) = 0;
W(N/2 + 2 : N)  = conj(W(N/2:-1:2));
w = ifft(W,'symmetric');
% figure;
% plot(time,w,'LineWidth',1.2);
% xlabel('Temps [sec]');
% ylabel('Composante fluctuante du vent w [m/s]');
% title('Processus al�atoire w(t)');
% 
%B) Impl�mentation d�un sch�ma d�int�gration pour r�soudre l��quation du mouvement et calcul de la r�ponse par ODE45:
F_b_omega = 1/4*rho*U(j)*B*[4*pi*W.', pi*B*W.']; % Force de turbulence en fr�quentiel
f_b_omega2 = ifft(F_b_omega,'symmetric'); % Devrait �tre �gal � f_b_omega car op�rations lin�aires => OK!
f_b_omega = 1/4*rho*U(j)*B*[4*pi*w.', pi*B*w.'];% Force de turbulence en temporel
f_b_omega = transpose(f_b_omega); %On transpose f_b_omega pour pouvoir extraire L_b et M_b comme �tant des vecteurs lignes
L_b = f_b_omega(1,:);
M_b = f_b_omega(2,:);
sigma_w = std(w); % Ecart-type de w(t)
sigma_Sw = sqrt(trapz(omega,Mat_S)*2); % Ecart_type de la PSD en fr�quentiel

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
VALEURS_PROPRES(:,j) = eig(M_final\A_final); % Pour la v�rification que les parties r�elles des valeurs propres soient n�gatives

f_b = zeros(10,length(omega));
f_b([3 4],:) = [L_b;M_b]; % Forces de turbulence en temporelle 

var = 1; % On d�marre � var = 2 pour �viter de calculer la fct de transfert quand omega = 0 (sinon on a une matrice remplie de valeurs infinies) 
while var <= length(omega)
H(:,:,var) = inv(1i*omega(var)*M_final - A_final); % FCT de transfert
H_zz(var) = H(1,1,var);
H_tt(var) = H(2,2,var);
var = var + 1;
end


% figure;
% plot(omega,abs(H_zz),'LineWidth',1.5); % H_zz
% xlim([0 2]);
% title('Fonction de r�ponse fr�quentielle H_{zz}(\omega)')
% xlabel('Pulsation \omega [rad/s]')
% ylabel('FRF')
% figure;
% plot(omega,abs(H_tt),'LineWidth',1.5); % H_thetatheta
% xlim([0 2]);
% title('Fonction de r�ponse fr�quentielle H_{\theta\theta}(\omega)')
% xlabel('Pulsation \omega [rad/s]')
% ylabel('FRF')
% 
tspan = time;
z0 = zeros(10,1); % structure initialement immobile
[t,z] = ode45(@(t,z) differentielle(t,z,A_final,M_final,f_b,time) ,tspan,z0); %ODE45 r�sout notre �quation diff�rentielle " M*z_dot = A*z + f_b"
% figure;
% plot(t,z(:,1),t,z(:,2),'LineWidth',1.5); % Graphe de la r�ponse z(t) et theta(t) de l'�quation du mouvement
% xlabel('Temps [sec]');
% title("R�ponse de l'�quation du mouvement par ODE45");
% legend('z(t)','\theta(t)');

%C) Calcul de la variance de la r�ponse:
dt = time(2)-time(1);

[PSD_z,F_z] = pwelch(z(:,1),[],[],length(omega),1/dt); % PSD de z(t)
F_zz = 2*pi*F_z; % Transformations de temporel vers fr�quentiel pour la comparaison dans le domaine fr�quentiel
PSD_zz = PSD_z/(4*pi);

[PSD_t,F_t] = pwelch(z(:,2),[],[],length(omega),1/dt); % PSD de theta(t)
F_tt = 2*pi*F_t; % Transformations de temporel vers fr�quentiel pour la comparaison dans le domaine fr�quentiel
PSD_tt = PSD_t/(4*pi);

% figure;
% semilogy(F_zz,PSD_zz,'LineWidth',1.5);
% xlim([0 2]);
% xlabel('Pulsation \omega [rad/s]');
% ylabel('PSD');
% title('Densit� spectrale de puissance S_{zz}(\omega)');
% figure;
% semilogy(F_tt,PSD_tt,'LineWidth',1.5);
% xlim([0 2]);
% xlabel('Pulsation \omega [rad/s]');
% ylabel('PSD');
% title('Densit� spectrale de puissance S_{\theta\theta}(\omega)');

% figure;
% semilogy(F_zz,PSD_zz,omega,Mat_S,'LineWidth',1.5); % Comparaison entre l'approche temporelle et l'approche fr�quentielle
% xlim([0 2]);
% xlabel('Pulsation \omega [rad/s]');
% ylabel('PSD');
% title('Densit� spectrale de puissance S_w(\omega)');
% legend('Approche temporelle','Approche fr�quentielle')

sigma_z(j) = std(z(:,1)); % Ecart-type de z(t)
sigma_theta(j) = std(z(:,2)); % Ecart-type de theta(t)

j = j + 1;
end

%D) Calcul de la vitesse de vent critique au flottement:

% U = 10:1:100;
% %U_crit = 76.86; % Valeur calcul�e dans la partie 1
% sigma_z(68:length(U)) = 0; % Etant donn� que le syst�me devient instable pass� la vitesse critique moyenne du vent U_crit, on va d�finir � z�ro les valeurs d'�cart-type pour U > U_crit
% sigma_theta(68:length(U)) = 0;
figure;
plot(U,sigma_z,'-o',U,sigma_theta*B/2,'-o','LineWidth',1.5);
xlabel('Vitesse moyenne horizontale U [m/s]');
ylabel('Ecart-type');
title('Ecart-type de la r�ponse');
legend('\sigma_{zz}','\sigma_{\theta\theta} * B/2')
