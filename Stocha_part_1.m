%% Partie 1 => Analyse fréquentielle du problème

function [omega_1,H_zz,H_tt,S_zz,S_tt] = Stocha_part_1(isComparison,U,Valeurs_utilisateur,Position)

if isComparison == 0
    TITRES = {'Pulsation(s) étudiée(s) \omega en [rad/s]: ','Vitesse(s) moyenne(s) horizontale(s) de vent étudiée(s) U en [m/s]: ','Vitesse(s) U à afficher dans [1 100]: (Si rien à afficher, [0])'};
    dlgtitle = 'Inputs pour analyse FRÉQUENTIELLE';
    dims = [1 100];
    Message_initial = {'0:0.01:2','10:1:100','[10]'};
    opts.Interpreter = 'tex';
    answer = inputdlg(TITRES,dlgtitle,dims,Message_initial,opts);
    omega_1 = str2num(answer{1});
    U = str2num(answer{2});
    Valeurs_utilisateur = str2num(answer{3});
    fprintf('Une des valeurs sélectionnées est U = %d [m/s]\n',Valeurs_utilisateur);
    [~,Position] = intersect(U,Valeurs_utilisateur);
    
elseif isComparison == 1
     TITRES = {'Pulsation(s) étudiée(s) \omega en [rad/s]: '};
    dlgtitle = 'Inputs pour analyse FRÉQUENTIELLE';
    dims = [1 100];
    Message_initial = {'0:0.01:2'};
    opts.Interpreter = 'tex';
    answer = inputdlg(TITRES,dlgtitle,dims,Message_initial,opts);
    omega_1 = str2num(answer{1});
end

if omega_1(1) < 1e-03
    omega_1(1) = 1e-03;
end
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
%omega_1 = 0:0.01:2;

%Propriétés du VENT incident:
%Masse volumique de l’air [kg/m^3]
rho = 1.22;
%Échelle de turbulence [m]
L_w = 20;
%Vitesse moyenne horizontale du vent [m/s]
%U = 10:1:100; % U = [0,100]: ATTENTION => Si vous générez une boucle de U, faites attention aux nombreuses figures qui seront plottées!!!
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

C_w =zeros(length(U),length(omega_1));
H = zeros(2,2,length(omega_1));
H_zz = zeros(length(U),length(omega_1)); H_tt = zeros(length(U),length(omega_1));
S_zz = zeros(length(U),length(omega_1)); S_tt = zeros(length(U),length(omega_1)); S_w = zeros(length(U),length(omega_1));
sigma_zz = zeros; sigma_tt = zeros; sigma_w = zeros;

j = 1;
while j <= length(U)
    
    i = 1;
    while i <= length(omega_1)
        [C_w(j,i),H(:,:,i),H_zz(j,i),H_tt(j,i)] = FRF(omega_1(i),U(j),B,rho,M,K,C);
        
        [S_zz(j,i),S_tt(j,i),S_w(j,i)] = PSD_omega_1(omega_1(i),U(j),B,rho,H(:,:,i),L_w,I_w);
        
        i = i + 1;
    end
    sigma_zz(j) = sqrt(trapz(omega_1,S_zz(j,:))*2); % Ecart-type de la PSD_zz 
    sigma_tt(j) = sqrt(trapz(omega_1,S_tt(j,:))*2); % Ecart-type de la PSD_thetatheta
    sigma_w(j) = sqrt(trapz(omega_1,S_w(j,:))*2);
        
    j = j + 1;
end

if isComparison == 0 
    var = 1;
    while var <= length(Valeurs_utilisateur)
        % Vérification Theodorsen:
        figure;
        hold on
        plot(omega_1,real(C_w(Position(var),:)),'LineWidth',1.5);
        plot(omega_1,imag(C_w(Position(var),:)),'LineWidth',1.5);
        hold off
        title('Fonction circulatoire de Theodorsen')
        xlabel('Pulsation \omega [rad/s]')
        legend('Partie réelle F(\omega)','Partie imaginaire G(\omega)');
        grid
        figure;
        plot(omega_1,abs(H_zz(Position(var),:)),'LineWidth',1.5); % Graphe H_zz
        title('Fonction de réponse fréquentielle H_{zz}(\omega)')
        xlabel('Pulsation \omega [rad/s]')
        ylabel('FRF')
        grid
        figure;
        plot(omega_1,abs(H_tt(Position(var),:)),'LineWidth',1.5); % Graphe H_\theta\theta
        title('Fonction de réponse fréquentielle H_{\theta\theta}(\omega)')
        xlabel('Pulsation \omega [rad/s]')
        ylabel('FRF')
        grid
        figure;
        semilogy(omega_1,abs(S_zz(Position(var),:)),'LineWidth',1.5); % Graphe PSD_zz
        title('Densité spectrale de puissance S_{zz}(\omega)')
        xlabel('Pulsation \omega [rad/s]')
        ylabel('PSD')
        grid
        figure;
        semilogy(omega_1,abs(S_tt(Position(var),:)),'LineWidth',1.5); % Graphe PSD_\theta\theta
        title('Densité spectrale de puissance S_{\theta\theta}(\omega)')
        xlabel('Pulsation \omega [rad/s]')
        ylabel('PSD')
        grid

        var = var + 1;
    end
end

if length(U) > 1
    U_crit = critical_speed(sigma_tt,sigma_zz,U);
    fprintf('La vitesse critique U_crit = %f [m/s]\n',U_crit);

    figure;
    hold on;
    plot(U,abs(sigma_zz),'-o','LineWidth',1.5);  % Graphe de l'écart-type sigma_zz
    plot(U,abs(sigma_tt*B/2),'-o','LineWidth',1.5); % Graphe de l'écart-type sigma_thetatheta
    xline(U_crit,'r','LineWidth',1.5);
    title('Ecart-type de la réponse');
    xlabel('Vitesse moyenne horizontale U [m/s]');
    ylabel('Ecart-type');
    legend('\sigma_{zz}','\sigma_{\theta\theta} * B/2','U_{critique}')
    grid
    hold off;
end

end

