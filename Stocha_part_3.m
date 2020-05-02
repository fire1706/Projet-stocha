function [] = Stocha_part_3()

TITRES = {'Vitesse(s) moyenne(s) horizontale(s) de vent étudiée(s) U en [m/s]: ','Voulez vous comparer approche FRÉQUENTIELLE et approche TEMPORELLE?: OUI => 1; NON => 0','Si OUI, pour quelle(s) valeur(s) de U en [m/s] ?: '};
dlgtitle = 'Inputs pour COMPARAISON';
dims = [1 100];
Message_initial = {'10:1:100','0 ou 1','10'};
answer = inputdlg(TITRES,dlgtitle,dims,Message_initial);
U = str2num(answer{1});
isComparison = str2num(answer{2});
if isComparison == 1
    fprintf('Vous avez décidé de comparer les deux approches\n');
elseif isComparison == 0
    fprintf('Vous avez décidé de ne pas comparer les deux approches\n');
end
U_plot = str2num(answer{3});
fprintf('Une des valeurs sélectionnées de comparaison est U = %d [m/s]\n',U_plot);
[~,Position] = intersect(U,U_plot);


[omega_1,H_zz,H_tt,S_zz,S_tt] = Stocha_part_1(isComparison,U,U_plot,Position);
[omega_2,H_zz2,H_tt2,F_zz,F_tt,PSD_zz,PSD_tt] = Stocha_part_2(isComparison,U,U_plot,Position);
if isComparison == 1
    
    figure;
    plot(omega_1,abs(H_zz),omega_2,abs(H_zz2),'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des FRF selon z')
    xlabel('Pulsation \omega [rad/s]')
    ylabel('FRF')
    legend('Approche fréquentielle','Approche temporelle')
    grid
    
    figure;
    plot(omega_1,abs(H_tt),omega_2,abs(H_tt2),'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des FRF selon \theta')
    xlabel('Pulsation \omega [rad/s]')
    ylabel('FRF')
    legend('Approche fréquentielle','Approche temporelle')
    grid
    
    figure;
    semilogy(omega_1,abs(S_zz),F_zz,PSD_zz,'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des PSD selon z')
    xlabel('Pulsation \omega [rad/s]')
    ylabel('PSD')
    legend('Approche fréquentielle','Approche temporelle')
    grid

    figure;
    semilogy(omega_1,abs(S_tt),F_tt,PSD_tt,'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des PSD selon \theta')
    xlabel('Pulsation \omega [rad/s]')
    ylabel('PSD')
    legend('Approche fréquentielle','Approche temporelle')
    grid
end


end