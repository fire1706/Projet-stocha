%% Partie 3 => Comparaison des donn�es fr�quentielles et temporelles

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

%Pulsation [rad/s]
w = 0:0.01:2;

%Propri�t�s du VENT incident:
%Masse volumique de l�air [kg/m^3]
rho = 1.22;
%�chelle de turbulence [m]
L_w = 20;
%Vitesse moyenne horizontale du vent [m/s]
U = 10:20:70; % U = [0,100]: ATTENTION => Si vous g�n�rez une boucle de U, faites attention aux nombreuses figures qui seront plott�es!!!
%Intensit� de turbulence [-]
I_w = 0.05;


MAX_sigma_tt = 0; % Variable pour la d�termination de U_critique (point D)
%% Partie 1

%A) Fct de transfert:
%Pulsation propre du tablier de pont [rad/s]
w_0 = [2*pi*f_z, 2*pi*f_z; 2*pi*f_theta, 2*pi*f_theta];
%Matrice de masse
M = [m,0;0,J];
%Matrice de raideur
K = w_0.^2.*M;
%Matrice d'amortissement
C = 2*xhi*sqrt(K*M); % Egalement valable: C = 2*xhi*w_0.*M;

j = 1;
while j <= length(U)
    
    H = zeros(2,2); % Initialisation fct de transfert
    i = 1;
    while i <= length(w)
        C_w(i) = C_omega(w(i),U(j),B); % Fct circulatoire de Theodorsen approxim�e par Jones
        
        q = pi*rho*U(j)^2*B;
        A = [B*w(i)^2/(4*U(j)^2) - C_w(i)*1i*w(i)/U(j) , B/(4*U(j))*1i*w(i) + C_w(i)*(1+B/(4*U(j))*1i*w(i)) ; -B/4*C_w(i)*1i*w(i)/U(j) , B/4*(B^2*w(i)^2/(32*U(j)^2) - B/(4*U(j))*1i*w(i) +  C_w(i)*(1+B/(4*U(j))*1i*w(i)))];
        
        H(:,:,i) = inv(-w(i)^2*M + 1i*w(i)*C + K - q*A); % Fct de transfert
        % H = X/fb c'est pour cela que fb n'est pas pris en compte dans notre fct de transfert;
        H_zz(i) = H(1,1,i);
        H_tt(i) = H(2,2,i);
        
        %B) Calcul densit� de puissance de la r�ponse:
        S_w(i) = S_omega(w(i),U(j)); % Calcul de la PSD de Von Karman
        F_b_omega = (1/4)*rho*U(j)*B*[4*pi;pi*B];
        S_b(:,:,i) = F_b_omega*S_w(i)*conj(transpose(F_b_omega));
        S_x(:,:,i) = H(:,:,i)*S_b(:,:,i)*conj(transpose(H(:,:,i)));
        S_zz(i) = S_x(1,1,i); % PSD_zz
        S_tt(i) = S_x(2,2,i); % PSD_thetatheta
        
        i = i + 1;
    end
    
    j = j + 1;
end

w1 = w;
%% Partie 2

%D�finition des variables pour cette partie:
%Temps maximal d�sir� [sec]
t_max = 600;
%Nombre d'�chantillons
N = 6000; % D�but: N = 2^(12) 
%Vecteur temps [sec]
time = linspace(0,t_max,N);

%Pas fr�quentiel [rad/s]
d_omega = 2*pi/t_max;
%Pulsation maximale d�sir�e [rad/s]
omega_max = (N-1)*d_omega;
%Pulsation [rad/s]
omega = (0:N-1)*d_omega;

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

%B) Impl�mentation d�un sch�ma d�int�gration pour r�soudre l��quation du mouvement et calcul de la r�ponse par ODE45:
f_b_omega = 1/4*rho*U(j)*B*[4*pi*w.', pi*B*w.'];% Force de turbulence en temporel
f_b_omega = transpose(f_b_omega);
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

p = 1;
while p <= length(time)
f_b(:,p) = [0;0;L_b(p);M_b(p);0;0;0;0;0;0]; % Forces de turbulence en temporelle 
p = p + 1;
end

var = 2; % On d�marre � var = 2 pour �viter de calculer la fct de transfert quand omega = 0 (sinon on a une matrice remplie de valeurs infinies) 
while var <= length(omega)
H2(:,:,var-1) = inv(1i*omega(var)*M_final - A_final); % FCT de transfert
H_zz2(var-1) = H2(1,1,var-1);
H_tt2(var-1) = H2(2,2,var-1);
var = var + 1;
end
valeurs_propres = real(eig(A_final)); % Pour la v�rification que les parties r�elles des valeurs propres soient n�gatives

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

tspan = time; % [0 t_max]
z0 = zeros(10,1); % structure initialement immobile
[t,z] = ode45(@(t,z) differentielle(t,z,A_final,M_final,f_b,time) ,tspan,z0); %ODE45 r�sout notre �quation diff�rentielle " M*z_dot = A*z + f_b"

%C) Calcul de la variance de la r�ponse:
dt = time(2)-time(1);

[PSD_Sw(:,j),F_Sw(:,j)] = pwelch(w,[],[],length(omega),1/dt);
F_Sw(:,j) = 2*pi*F_Sw(:,j);
PSD_Sw(:,j) = PSD_Sw(:,j)/(4*pi);

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

j = j + 1;

end

%% G�n�ration des graphiques pour les comparaisons:

figure;
semilogy(F_zz,PSD_zz,w1,abs(S_zz),'LineWidth',1.5)
xlim([0 2]);
title('Comparaison des PSD selon z')
xlabel('Pulsation \omega [rad/s]')
ylabel('PSD')
legend('Approche temporelle','Approche fr�quentielle')

figure;
semilogy(F_tt,PSD_tt,w1,abs(S_tt),'LineWidth',1.5)
xlim([0 2]);
title('Comparaison des PSD selon \theta')
xlabel('Pulsation \omega [rad/s]')
ylabel('PSD')
legend('Approche temporelle','Approche fr�quentielle')


figure;
var = 1;
while var <= length(U)
hold on
plot(F_Sw(:,var),log10(PSD_Sw(:,var)),'LineWidth',1.5); % Comparaison entre l'approche temporelle et l'approche fr�quentielle
plot(omega,log10(Mat_S),'LineWidth',1.5);
xlim([0 2]);
xlabel('Pulsation \omega [rad/s]');
ylabel('PSD');
title('Densit� spectrale de puissance S_w(\omega)');
legend('Approche temporelle','Approche fr�quentielle');
hold off
var = var + 1;
end

% figure;
% plot(omega,H_zz2,w1,abs(H_zz),'LineWidth',1.5)
% xlim([0 2]);
% title('Comparaison des FRF selon z')
% xlabel('Pulsation \omega [rad/s]')
% ylabel('FRF')
% legend('Approche fr�quentielle','Approche temporelle')
% 
% figure;
% plot(omega,H_tt2,w1,abs(H_tt),'LineWidth',1.5)
% xlim([0 2]);
% title('Comparaison des FRF selon \theta')
% xlabel('Pulsation \omega [rad/s]')
% ylabel('FRF')
% legend('Approche fr�quentielle','Approche temporelle')
