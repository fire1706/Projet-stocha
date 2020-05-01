function [sigma_zz,sigma_tt,sigma_w,sigma_z,sigma_theta,sigma_Sw] = Stocha_part_3(omega_1,t_max,N,U,U_plot,Position,Flow,isComparison2)

% TITRES = {'Écoulement est instationnaire ou quasi-stationnaire?:  Si insta => 0; si quasi => 1','Vitesse(s) moyenne(s) horizontale(s) de vent étudiée(s) U en [m/s]: ','Voulez vous comparer approche FRÉQUENTIELLE et approche TEMPORELLE?: OUI => 1; NON => 0','Si OUI, pour quelle(s) valeur(s) de U en [m/s] ?: '};
% dlgtitle = 'Inputs';
% dims = [1 100];
% Message_initial = {'0 ou 1','10:1:100','0 ou 1','10'};
% answer = inputdlg(TITRES,dlgtitle,dims,Message_initial);
% Flow = str2num(answer{1});
% U = str2num(answer{2});
% isComparison2 = str2num(answer{3});
% if isComparison2 == 1
%     fprintf('Vous avez décidé de comparer les deux approches\n');
% elseif isComparison2 == 0
%     fprintf('Vous avez décidé de ne pas comparer les deux approches\n');
% end
% U_plot = str2num(answer{4});
% fprintf('Une des valeurs sélectionnées de comparaison est U = %d [m/s]\n',U_plot);
% [~,Position] = intersect(U,U_plot);


[omega_1,H_zz,H_tt,S_zz,S_tt,sigma_zz,sigma_tt,sigma_w] = Stocha_part_1(omega_1,U,U_plot,Position,Flow,isComparison2);
[omega_2,H_zz2,H_tt2,F_zz,F_tt,F_ww,PSD_zz,PSD_tt,PSD_ww,S_w,sigma_z,sigma_theta,sigma_Sw] = Stocha_part_2(t_max,N,U,U_plot,Position,Flow,isComparison2);
if isComparison2 == 1
    loop = 1;
    while loop <= length(U_plot)
    figure;
    plot(omega_1,abs(H_zz(Position(loop),:)),omega_2,abs(H_zz2(Position(loop),:)),'LineWidth',1.5)
    xlim([0 2]);
    title(['Comparaison des FRF selon z pour U = ',num2str(U_plot(loop)),' [mm/s]']);
    xlabel('Pulsation \omega [rad/s]')
    ylabel('FRF')
    legend('Approche fréquentielle','Approche temporelle')
    grid
    
    figure;
    plot(omega_1,abs(H_tt(Position(loop),:)),omega_2,abs(H_tt2(Position(loop),:)),'LineWidth',1.5)
    xlim([0 2]);
    title(['Comparaison des FRF selon \theta pour U = ',num2str(U_plot(loop)),' [mm/s]']);
    xlabel('Pulsation \omega [rad/s]')
    ylabel('FRF')
    legend('Approche fréquentielle','Approche temporelle')
    grid
    
    figure;
    semilogy(omega_1,abs(S_zz(Position(loop),:)),F_zz(:,Position(loop)),PSD_zz(:,Position(loop)),'LineWidth',1.5)
    xlim([0 2]);
    title(['Comparaison des PSD selon z pour U = ',num2str(U_plot(loop)),' [mm/s]']);
    xlabel('Pulsation \omega [rad/s]')
    ylabel('PSD')
    legend('Approche fréquentielle','Approche temporelle')
    grid

    figure;
    semilogy(omega_1,abs(S_tt(Position(loop),:)),F_tt(:,Position(loop)),PSD_tt(:,Position(loop)),'LineWidth',1.5)
    xlim([0 2]);
    title(['Comparaison des PSD selon \theta pour U = ',num2str(U_plot(loop)),' [mm/s]']);
    xlabel('Pulsation \omega [rad/s]')
    ylabel('PSD')
    legend('Approche fréquentielle','Approche temporelle')
    grid
    
    figure;
    semilogy(omega_2,S_w(Position(loop),:),F_ww(:,Position(loop)),PSD_ww(:,Position(loop)),'LineWidth',1.5)
    xlim([0 2]);
    title(['Comparaison des densités spectrales de puissance S_w(\omega) pour U = ',num2str(U_plot(loop)),' [mm/s]']);
    xlabel('Pulsation \omega [rad/s]')
    ylabel('PSD')
    legend('Approche fréquentielle','Approche temporelle')
    grid
 
    loop = loop + 1;
    end
end


end