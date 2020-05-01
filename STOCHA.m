function [sigma_zz,sigma_tt,sigma_w,sigma_z,sigma_theta,sigma_Sw] = STOCHA()

%% Quasi-Statio/Instatio:

TITRES = {'Voulez-vous comparer un écoulement quasi-stationnaire avec un instationnaire?: Si OUI => 1  Si NON => 0','Si NON, écoulement est instationnaire ou quasi-stationnaire?:  Si quasi => 1  Si insta => 0'};
dlgtitle = 'Inputs';
dims = [1 150];
Message_initial = {'0 ou 1','0 ou 1'};
answer = inputdlg(TITRES,dlgtitle,dims,Message_initial);
isComparison1 = str2num(answer{1});
Flow = str2num(answer{2});
if isComparison1 == 1
    fprintf('Vous avez décidé de comparer les deux écoulements\n');
elseif isComparison1 == 0
    fprintf('Vous avez décidé de ne pas comparer les deux écoulements\n');
end

%% Temp/Freq

if isComparison1 == 0
    
    TITRES = {'Voulez vous comparer approche FRÉQUENTIELLE et approche TEMPORELLE?:  Si OUI => 1  Si NON => 0','Si NON, approche FREQUENTIELLE ou TEMPORELLE?:  Si freq => 1  Si temp => 0'};
    dlgtitle = 'Inputs';
    dims = [1 150];
    Message_initial = {'0 ou 1','0 ou 1'};
    answer = inputdlg(TITRES,dlgtitle,dims,Message_initial);
    isComparison2 = str2num(answer{1});
    Approach = str2num(answer{2});
    if isComparison2 == 1
        fprintf('Vous avez décidé de comparer les deux approches\n');
    elseif isComparison2 == 0
        fprintf('Vous avez décidé de ne pas comparer les deux approches\n');
    end
end

%% CODES 

if isComparison1 == 1
    
    TITRES = {'Pulsation(s) étudiée(s) \omega en [rad/s] pour approche FRÉQUENTIELLE: ','Temps maximal t_{max} étudié en [sec] pour approche TEMPORELLE: ','Nombre échantillons N générés [-] pour approche TEMPORELLE: ','Vitesse(s) moyenne(s) horizontale(s) de vent étudiée(s) U en [m/s]: ','Vitesse(s) U à afficher dans [1 100]: '};
    dlgtitle = 'Inputs pour analyse FRÉQUENTIELLE';
    dims = [1 150];
    Message_initial = {'0:0.01:2','600','600','10:1:100','[10]'};
    opts.Interpreter = 'tex';
    answer = inputdlg(TITRES,dlgtitle,dims,Message_initial,opts);
    omega_1 = str2num(answer{1});
    t_max = str2num(answer{2});
    N = str2num(answer{3});
    U = str2num(answer{4});
    Valeurs_utilisateur = str2num(answer{5});
    fprintf('Une des valeurs sélectionnées est U = %d [m/s]\n',Valeurs_utilisateur);
    [~,Position] = intersect(U,Valeurs_utilisateur);
    isComparison2 = 1;
    i = 1;
    while i <= 2
        if i == 1
            Flow = 0;
            [omega_1(:,:,i),H_zz(:,:,i),H_tt(:,:,i),S_zz(:,:,i),S_tt(:,:,i)] = Stocha_part_1(omega_1,U,Valeurs_utilisateur,Position,Flow,isComparison2);
            [omega_2(:,:,i),H_zz2(:,:,i),H_tt2(:,:,i),F_zz(:,:,i),F_tt(:,:,i),F_ww(:,:,i),PSD_zz(:,:,i),PSD_tt(:,:,i),PSD_ww(:,:,i),S_w(:,:,i)] = Stocha_part_2(t_max,N,U,Valeurs_utilisateur,Position,Flow,isComparison2);
        elseif i == 2
            Flow = 1;
            [omega_1(:,:,i),H_zz(:,:,i),H_tt(:,:,i),S_zz(:,:,i),S_tt(:,:,i)] = Stocha_part_1(omega_1,U,Valeurs_utilisateur,Position,Flow,isComparison2);
            [omega_2(:,:,i),H_zz2(:,:,i),H_tt2(:,:,i),F_zz(:,:,i),F_tt(:,:,i),F_ww(:,:,i),PSD_zz(:,:,i),PSD_tt(:,:,i),PSD_ww(:,:,i),S_w(:,:,i)] = Stocha_part_2(t_max,N,U,Valeurs_utilisateur,Position,Flow,isComparison2);
        end
        i = i + 1;
    end
elseif isComparison1 == 0
    
    if isComparison2 == 0
        if Approach == 0
            
            TITRES = {'Temps maximal t_{max} étudié en [sec] pour approche TEMPORELLE: ','Nombre échantillons N générés [-] pour approche TEMPORELLE: ','Vitesse(s) moyenne(s) horizontale(s) de vent étudiée(s) U en [m/s]: ','Vitesse(s) U à afficher dans [1 100]: '};
            dlgtitle = 'Inputs pour analyse FRÉQUENTIELLE';
            dims = [1 150];
            Message_initial = {'600','600','10:1:100','[10]'};
            opts.Interpreter = 'tex';
            answer = inputdlg(TITRES,dlgtitle,dims,Message_initial,opts);
            t_max = str2num(answer{1});
            N = str2num(answer{2});
            U = str2num(answer{3});
            Valeurs_utilisateur = str2num(answer{4});
            fprintf('Une des valeurs sélectionnées est U = %d [m/s]\n',Valeurs_utilisateur);
            [~,Position] = intersect(U,Valeurs_utilisateur);
            
            Stocha_part_2(t_max,N,U,Valeurs_utilisateur,Position,Flow,isComparison2);
            
        elseif Approach == 1
            
            TITRES = {'Pulsation(s) étudiée(s) \omega en [rad/s]: ','Vitesse(s) moyenne(s) horizontale(s) de vent étudiée(s) U en [m/s]: ','Vitesse(s) U à afficher dans [1 100]: '};
            dlgtitle = 'Inputs pour analyse FRÉQUENTIELLE';
            dims = [1 150];
            Message_initial = {'0:0.01:2','10:1:100','[10]'};
            opts.Interpreter = 'tex';
            answer = inputdlg(TITRES,dlgtitle,dims,Message_initial,opts);
            omega_1 = str2num(answer{1});
            U = str2num(answer{2});
            Valeurs_utilisateur = str2num(answer{3});
            fprintf('Une des valeurs sélectionnées est U = %d [m/s]\n',Valeurs_utilisateur);
            [~,Position] = intersect(U,Valeurs_utilisateur);

            Stocha_part_1(omega_1,U,Valeurs_utilisateur,Position,Flow,isComparison2);
            
        end
    elseif isComparison2 == 1
        
        TITRES = {'Pulsation(s) étudiée(s) \omega en [rad/s] pour approche FRÉQUENTIELLE: ','Temps maximal t_{max} étudié en [sec] pour approche TEMPORELLE: ','Nombre échantillons N générés [-] pour approche TEMPORELLE: ','Vitesse(s) moyenne(s) horizontale(s) de vent étudiée(s) U en [m/s]: ','Vitesse(s) U à afficher dans [1 100]: '};
        dlgtitle = 'Inputs pour analyse FRÉQUENTIELLE';
        dims = [1 150];
        Message_initial = {'0:0.01:2','600','600','10:1:100','[10]'};
        opts.Interpreter = 'tex';
        answer = inputdlg(TITRES,dlgtitle,dims,Message_initial,opts);
        omega_1 = str2num(answer{1});
        t_max = str2num(answer{2});
        N = str2num(answer{3});
        U = str2num(answer{4});
        Valeurs_utilisateur = str2num(answer{5});
        fprintf('Une des valeurs sélectionnées est U = %d [m/s]\n',Valeurs_utilisateur);
        [~,Position] = intersect(U,Valeurs_utilisateur);
        
       [sigma_zz,sigma_tt,sigma_w,sigma_z,sigma_theta,sigma_Sw] = Stocha_part_3(omega_1,t_max,N,U,Valeurs_utilisateur,Position,Flow,isComparison2);
        
    end
end

%% PLOT si comparaison quasi et insta:
if isComparison1 == 1
    
    figure;
    plot(omega_1(:,:,1),abs(H_zz(:,:,1)),omega_1(:,:,2),abs(H_zz(:,:,2)),'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des FRF selon z par approche fréquentielle');
    xlabel('Pulsation \omega [rad/s]')
    ylabel('FRF')
    legend('Écoulement instationnaire','Écoulement quasi-stationnaire');
    grid

    figure;
    plot(omega_2(:,:,1),abs(H_zz2(:,:,1)),omega_2(:,:,2),abs(H_zz2(:,:,2)),'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des FRF selon z par approche temporelle');
    xlabel('Pulsation \omega [rad/s]')
    ylabel('FRF')
    legend('Écoulement instationnaire','Écoulement quasi-stationnaire');
    grid


    figure;
    plot(omega_1(:,:,1),abs(H_tt(:,:,1)),omega_1(:,:,2),abs(H_tt(:,:,2)),'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des FRF selon \theta par approche fréquentielle');
    xlabel('Pulsation \omega [rad/s]')
    ylabel('FRF')
    legend('Écoulement instationnaire','Écoulement quasi-stationnaire');
    grid

    figure;
    plot(omega_2(:,:,1),abs(H_tt2(:,:,1)),omega_2(:,:,2),abs(H_tt2(:,:,2)),'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des FRF selon \theta par approche temporelle');
    xlabel('Pulsation \omega [rad/s]')
    ylabel('FRF')
    legend('Écoulement instationnaire','Écoulement quasi-stationnaire');
    grid

    figure;
    semilogy(omega_1(:,:,1),abs(S_zz(:,:,1)),omega_1(:,:,2),abs(S_zz(:,:,2)),'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des PSD selon z par approche fréquentielle');
    xlabel('Pulsation \omega [rad/s]')
    ylabel('PSD')
    legend('Écoulement instationnaire','Écoulement quasi-stationnaire');
    grid

    figure;
    semilogy(F_zz(:,:,1),PSD_zz(:,:,1),F_zz(:,:,2),PSD_zz(:,:,2),'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des PSD selon z par approche temporelle');
    xlabel('Pulsation \omega [rad/s]')
    ylabel('PSD')
    legend('Écoulement instationnaire','Écoulement quasi-stationnaire');
    grid

    figure;
    semilogy(omega_1(:,:,1),abs(S_tt(:,:,1)),omega_1(:,:,2),abs(S_tt(:,:,2)),'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des PSD selon \theta par approche fréquentielle');
    xlabel('Pulsation \omega [rad/s]')
    ylabel('PSD')
    legend('Écoulement instationnaire','Écoulement quasi-stationnaire');
    grid

    figure;
    semilogy(F_tt(:,:,1),PSD_tt(:,:,1),F_tt(:,:,2),PSD_tt(:,:,2),'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des PSD selon \theta par approche temporelle');
    xlabel('Pulsation \omega [rad/s]')
    ylabel('PSD')
    legend('Écoulement instationnaire','Écoulement quasi-stationnaire');
    grid


    figure;
    semilogy(omega_2(:,:,1),S_w(:,:,1),F_ww(:,:,1),PSD_ww(:,:,1),'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des PSD S_w(\omega) pour écoulement instationnaire');
    xlabel('Pulsation \omega [rad/s]')
    ylabel('PSD')
    legend('Approche fréquentielle','Approche temporelle');
    grid

    figure;
    semilogy(omega_2(:,:,2),S_w(:,:,2),F_ww(:,:,2),PSD_ww(:,:,2),'LineWidth',1.5)
    xlim([0 2]);
    title('Comparaison des PSD S_w(\omega) pour écoulement quasi-stationnaire');
    xlabel('Pulsation \omega [rad/s]')
    ylabel('PSD')
    legend('Approche fréquentielle','Approche temporelle');
    grid
    
end

end
       
    
