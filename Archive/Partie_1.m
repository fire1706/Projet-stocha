%% Partie 1 => Analyse fréquentielle du problème
%clear all;
%% DONNEES
%Propriétés structurelles du PONT:
%Masse par unité de longueur [kg/m]
m = 22740;
%Moment d’inertie par unité de longueur [kg.m^2/m]
J = 2.47*10^6;
%Largeur du tablier [m]
B = 31;
%Amortissement structurel [-]
xhi = 0.003;
%Fréquence propre verticale [Hz]
f_z = 0.10;
%Fréquence propre de torsion [Hz]
f_theta = 0.278;

%Pulsation [rad/s]
w = 0:0.01:2;

%Propriétés du VENT incident:
%Masse volumique de l’air [kg/m^3]
rho = 1.22;
%Échelle de turbulence [m]
L_w = 20;
%Vitesse moyenne horizontale du vent [m/s]
U = 10:1:100; % U = [0,100]: ATTENTION => Si vous générez une boucle de U, faites attention aux nombreuses figures qui seront plottées!!!
%Intensité de turbulence [-]
I_w = 0.05;

%% Résolution

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
        C_w(i) = C_omega(w(i),U(j),B); % Fct circulatoire de Theodorsen approximée par Jones
        
        q = pi*rho*U(j)^2*B;
        A = [B*w(i)^2/(4*U(j)^2) - C_w(i)*1i*w(i)/U(j) , B/(4*U(j))*1i*w(i) + C_w(i)*(1+B/(4*U(j))*1i*w(i)) ; -B/4*C_w(i)*1i*w(i)/U(j) , B/4*(B^2*w(i)^2/(32*U(j)^2) - B/(4*U(j))*1i*w(i) +  C_w(i)*(1+B/(4*U(j))*1i*w(i)))];
        
        H(:,:,i) = inv(-w(i)^2*M + 1i*w(i)*C + K - q*A); % Fct de transfert
        % H = X/fb c'est pour cela que fb n'est pas pris en compte dans notre fct de transfert;
        H_zz(i) = H(1,1,i);
        H_tt(i) = H(2,2,i);
        
        %B) Calcul densité de puissance de la réponse:
        S_w(i) = S_omega(w(i),U(j)); % Calcul de la PSD de Von Karman
        F_b_omega = (1/4)*rho*U(j)*B*[4*pi;pi*B];
        S_b(:,:,i) = F_b_omega*S_w(i)*conj(transpose(F_b_omega));
        S_x(:,:,i) = H(:,:,i)*S_b(:,:,i)*conj(transpose(H(:,:,i)));
        S_zz(i) = S_x(1,1,i); % PSD_zz
        S_tt(i) = S_x(2,2,i); % PSD_thetatheta
        
        i = i + 1;
    end
    
% figure;
% grid;
% plot(w,abs(H_zz),'LineWidth',1.5); % Graphe H_zz
% title('Fonction de réponse fréquentielle H_{zz}(\omega)')
% xlabel('Pulsation \omega [rad/s]')
% ylabel('FRF')
% figure;
% grid;
% plot(w,abs(H_tt),'LineWidth',1.5); % Graphe H_thetatheta
% title('Fonction de réponse fréquentielle H_{\theta\theta}(\omega)')
% xlabel('Pulsation \omega [rad/s]')
% ylabel('FRF')
% 
% figure;
% grid;
% semilogy(w,abs(S_zz),'LineWidth',1.5); % Graphe PSD_zz
% title('Densité spectrale de puissance S_{zz}(\omega)')
% xlabel('Pulsation \omega [rad/s]')
% ylabel('PSD')
% figure;
% grid;
% semilogy(w,abs(S_tt),'LineWidth',1.5); % Graphe PSD_thetatheta
% title('Densité spectrale de puissance S_{\theta\theta}(\omega)')
% xlabel('Pulsation \omega [rad/s]')
% ylabel('PSD')

%C) Calcul de la variance de la réponse:
sigma_zz(j) = sqrt(trapz(w,S_zz)); % Ecart-type de la PSD_zz 
sigma_tt(j) = sqrt(trapz(w,S_tt)); % Ecart-type de la PSD_thetatheta
   
    % Vérification Theodorsen:
%     figure;
%     hold on
%     plot(w,real(C_w),'LineWidth',1.5);
%     plot(w,imag(C_w),'LineWidth',1.5);
%     hold off
%     title('Fonction circulatoire de Theodorsen')
%     xlabel('Pulsation \omega [rad/s]')
%     legend('Partie réelle F(\omega)','Partie imaginaire G(\omega)');
    
    j = j + 1;
end

%D) Calcul de la vitesse du vent critique au flottement:
U_crit = 0;
cumul_MAXtt = cummax(sigma_tt);
cumul_MAXzz = cummax(sigma_zz);
i = 2;
while i <= length(U)-1
    if cumul_MAXtt(i-1) == cumul_MAXtt(i) && cumul_MAXzz(i-1) == cumul_MAXzz(i)
        U_crit = U(i-1);
        break;
    end
    i = i + 1;
end

% figure;
% grid;
% hold on;
% plot(U,abs(sigma_zz),'-o','LineWidth',1.5);  % Graphe de l'écart-type sigma_zz
% plot(U,abs(sigma_tt*B/2),'-o','LineWidth',1.5); % Graphe de l'écart-type sigma_thetatheta
% xline(U_crit,'r','LineWidth',1.5);
% title('Ecart-type de la réponse');
% xlabel('Vitesse moyenne horizontale U [m/s]');
% ylabel('Ecart-type');
% legend('\sigma_{zz}','\sigma_{\theta\theta} * B/2','U_{critique}')
% hold off;
