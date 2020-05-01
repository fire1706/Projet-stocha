%% Partie 2 => Analyse temporelle du problème

function [omega_2,H_zz2,H_tt2,F_zz,F_tt,F_ww,PSD_zz,PSD_tt,PSD_ww,S_w,sigma_z,sigma_theta,sigma_Sw] = Stocha_part_2(t_max,N,U,Valeurs_utilisateur,Position,Flow,isComparison2)

% if isComparison2 == 0
%     TITRES = {'Temps maximal t_max étudié en [sec]: ','Nombre échantillons N générés [-]: ','Vitesse(s) moyenne(s) horizontale(s) de vent étudiée(s) U en [m/s]: ','Vitesse(s) U à afficher dans [1 100]: (Si rien à afficher, [0])'};
%     dlgtitle = 'Inputs pour analyse TEMPORELLE';
%     dims = [1 100];
%     Message_initial = {'600','600','10:1:70','[10]'};
%     answer = inputdlg(TITRES,dlgtitle,dims,Message_initial);
%     t_max = str2num(answer{1});
%     N = str2num(answer{2});
%     U = str2num(answer{3});
%     Valeurs_utilisateur = str2num(answer{4});
%     fprintf('Une des valeurs sélectionnées est U = %d [m/s]\n',Valeurs_utilisateur);
%     [~,Position] = intersect(U,Valeurs_utilisateur);
% 
% elseif isComparison2 == 1
%     TITRES = {'Temps maximal t_max étudié en [sec]: ','Nombre échantillons N générés [-]: '};
%     dlgtitle = 'Inputs pour analyse TEMPORELLE';
%     dims = [1 100];
%     Message_initial = {'600','600'};
%     answer = inputdlg(TITRES,dlgtitle,dims,Message_initial);
%     t_max = str2num(answer{1});
%     N = str2num(answer{2});
% end
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

%Propriétés du VENT incident:
%Masse volumique de l’air [kg/m^3]
rho = 1.22;
%Échelle de turbulence [m]
L_w = 20;

%Intensité de turbulence [-]
I_w = 0.05;

%Valeurs utiles trouvées à la Partie 1:
%Pulsation propre du tablier de pont [rad/s]
w_0 = [2*pi*f_z,2*pi*f_z;2*pi*f_theta,2*pi*f_theta];
%Matrice de masse
M = [m,0;0,J];
%Matrice de raideur
K = w_0.^2.*M;
%Matrice d'amortissement
C = 2*xhi*sqrt(K*M);% Egalement valable: C = 2*xhi*w_0.*M;

%Définition des variables pour cette partie:
%Temps maximal désiré [sec]
%t_max = 600;
%Nombre d'échantillons
%N = 600; % Début: N = 2^(12) 
%Vecteur temps [sec]
time = linspace(0,t_max,N);

%Pas fréquentiel [rad/s]
d_omega = 2*pi/t_max;
%Pulsation maximale désirée [rad/s]
omega_max = (N-1)*d_omega;
%Pulsation [rad/s]
omega_2 = (0:N-1)*d_omega;
omega_2(1) = 1e-03;

%% Résolution
w = zeros(length(U),length(omega_2));
Zout = zeros(length(omega_2),10,length(U)); VALEURS_PROPRES = zeros(10,length(U));
%PSD_zz = zeros; PSD_tt = zeros; 
H_zz2 = zeros(length(U),length(omega_2)); H_tt2 = zeros(length(U),length(omega_2));

j = 1;
while j <= length(U)
    [w(j,:),S_w(j,:),W(j,:)] = Random_process(d_omega,omega_2,N,U(j),L_w,I_w);
    [t(:,j),Zout(:,:,j),VALEURS_PROPRES(:,j),A_final,M_final] = Diff_equation_for_ode45(omega_2,time,U(j),B,rho,W,w(j,:),M,K,C,Flow);
    [PSD_zz(:,j),F_zz(:,j),PSD_tt(:,j),F_tt(:,j),PSD_ww(:,j),F_ww(:,j)] = PSD_omega_2(time,omega_2,Zout(:,:,j),w(j,:));
    i = 1;
    while i <= length(omega_2)
        [H_zz2(j,i),H_tt2(j,i)] = FRF_part2(omega_2(i),A_final,M_final);
        i = i + 1;
    end
    
    sigma_w(j) = std(w(j,:)); % Ecart-type de w(t)
    sigma_Sw(j) = sqrt(trapz(omega_2,S_w(j,:))*2); % Ecart_type de la PSD en fréquentiel
    sigma_z(j) = std(Zout(:,1,j)); % Ecart-type de z(t)
    sigma_theta(j) = std(Zout(:,2,j)); % Ecart-type de theta(t)
    j = j + 1;
    
end

if isComparison2 == 0
    loop = 1;
    while loop <= length(Valeurs_utilisateur)
        figure;
        plot(time,w(Position(loop),:),'LineWidth',1.2);
        xlabel('Temps [sec]');
        ylabel('Composante fluctuante du vent w [m/s]');
        title(['Processus aléatoire w(t) pour U = ',num2str(Valeurs_utilisateur(loop)),' [mm/s]']);
        grid
        figure;
        plot(t(:,Position(loop)),Zout(:,1,Position(loop)),t(:,Position(loop)),Zout(:,2,Position(loop)),'LineWidth',1.5); % Graphe de la réponse z(t) et theta(t) de l'équation du mouvement
        xlabel('Temps [sec]');
        title(["Réponse de l'équation du mouvement par ODE45 pour U = ",num2str(Valeurs_utilisateur(loop))," [mm/s]"]);
        legend('z(t)','\theta(t)');
        grid
        figure;
        plot(omega_2,abs(H_zz2(Position(loop),:)),'LineWidth',1.5); % Graphe H_zz
        title(['Fonction de réponse fréquentielle H_{zz}(\omega) pour U = ',num2str(Valeurs_utilisateur(loop)),' [mm/s]']);
        xlabel('Pulsation \omega [rad/s]')
        ylabel('FRF')
        xlim([0 2]);
        grid
        figure;
        plot(omega_2,abs(H_tt2(Position(loop),:)),'LineWidth',1.5); % Graphe H_\theta\theta
        title(['Fonction de réponse fréquentielle H_{\theta\theta}(\omega) pour U = ',num2str(Valeurs_utilisateur(loop)),' [mm/s]']);
        xlabel('Pulsation \omega [rad/s]')
        ylabel('FRF')
        xlim([0 2]);
        grid
        figure;
        semilogy(F_zz(:,Position(loop)),PSD_zz(:,Position(loop)),'LineWidth',1.5);
        xlim([0 2]);
        xlabel('Pulsation \omega [rad/s]');
        ylabel('PSD');
        title(['Densité spectrale de puissance S_{zz}(\omega) pour U = ',num2str(Valeurs_utilisateur(loop)),' [mm/s]']);
        grid
        figure;
        semilogy(F_tt(:,Position(loop)),PSD_tt(:,Position(loop)),'LineWidth',1.5);
        xlim([0 2]);
        xlabel('Pulsation \omega [rad/s]');
        ylabel('PSD');
        title(['Densité spectrale de puissance S_{\theta\theta}(\omega) pour U = ',num2str(Valeurs_utilisateur(loop)),' [mm/s]']);
        grid

        loop = loop + 1;
    end
end

if length(U) > 1
    [U_crit,sigma_z,sigma_theta] = critical_speed_part2(VALEURS_PROPRES,U,sigma_z,sigma_theta);
    figure;
    hold on
    plot(U,sigma_z,'-o',U,sigma_theta*B/2,'-o','LineWidth',1.5);
    xline(U_crit,'r','LineWidth',1.5);
    xlabel('Vitesse moyenne horizontale U [m/s]');
    ylabel('Ecart-type');
    title('Ecart-type de la réponse');
    legend('\sigma_{zz}','\sigma_{\theta\theta} * B/2','U critique')
    grid
    hold off
    fprintf('La vitesse critique U_crit est de %f [m/s]\n',U_crit);
end

end