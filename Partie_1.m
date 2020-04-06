%% Partie 1 => Analyse fréquentielle du problème
clear all;
%% DONNEES
%Propriétés du PONT:
i = 1;
w = 0:0.01:2;
m = 22740;
J = 2.47*10^6;
B = 31;
xhi = 0.003;
f_z = 0.10;
f_theta = 0.278;

%Propriétés du VENT:
rho = 1.22;
L_w = 20;
U = 10:1:100; %U = [0,100]: ATTENTION => Si vous générez une boucle de U, faites attention aux nombreuses figures qui seront plottées!!!
I_w = 0.05;


MAX_sigma_tt = 0; % Variable pour la détermination de U_critique (point D)
%% Résolution

%A) Fct de transfert:
w_0 = [2*pi*f_z, 2*pi*f_z; 2*pi*f_theta, 2*pi*f_theta];
M = [m,0;0,J];
K = w_0.^2.*M;
C = 2*xhi*sqrt(K*M); % Egalement valable: C = 2*xhi*w_0.*M;

j = 1;
while j <= length(U)
    
    % Approximation par Jones:
    a0 = 1;
    a1 = -0.165;
    a2 = -0.335;
    b1 = 0.0455;
    b2 = 0.3;
    w_1 = 2*U(j)*b1/B;
    w_2 = 2*U(j)*b2/B;
    
    H = zeros(2,2);
    i = 1;
    while i <= length(w)
        C_w(i) = C_omega(w(i),U(j),B); % Fct circulatoire de Theodorsen
        
        q = pi*rho*U(j)^2*B;
        A = [B*w(i)^2/(4*U(j)^2) - C_w(i)*1i*w(i)/U(j) , B/(4*U(j))*1i*w(i) + C_w(i)*(1+B/(4*U(j))*1i*w(i)) ; -B/4*C_w(i)*1i*w(i)/U(j) , B/4*(B^2*w(i)^2/(32*U(j)^2) - B/(4*U(j))*1i*w(i) +  C_w(i)*(1+B/(4*U(j))*1i*w(i)))];
        
        H(:,:,i) = inv(-w(i)^2*M + 1i*w(i)*C + K - q*A); % Fct de transfert
        % H = X/fb c'est pour cela que fb n'est pas pris en compte dans notre fct de transfert;
        H_zz(i) = H(1,1,i);
        H_tt(i) = H(2,2,i);
        
        %B) Calcul densité de puissance de la réponse:
        S_w(i) = S_omega(w(i),U(j));
        F_b_omega = (1/4)*rho*U(j)*B*[4*pi;pi*B];
        S_b(:,:,i) = F_b_omega*S_w(i)*conj(transpose(F_b_omega));
        S_x(:,:,i) = H(:,:,i)*S_b(:,:,i)*conj(transpose(H(:,:,i)));
        S_zz(i) = S_x(1,1,i); % PSD_zz
        S_tt(i) = S_x(2,2,i); % PSD_thetatheta
        
        i = i + 1;
    end
    
% figure;
% plot(w,abs(H_zz),'LineWidth',1.5); % Graphe H_zz
% title( 'Nodal FRF H_{zz}(\omega)' )
% xlabel('Pulsations \omega [rad/s]')
% ylabel('FRF')
% figure;
% plot(w,abs(H_tt),'LineWidth',1.5); % Graphe H_thetatheta
% title( 'Nodal FRF H_{\theta\theta}(\omega)' )
% xlabel('Pulsations \omega [rad/s]')
% ylabel('FRF')
% figure;
% 
% semilogy(w,abs(S_zz),'LineWidth',1.5); % Graphe PSD_zz
% title( 'Nodal PSD S_{zz}(\omega)' )
% xlabel('Pulsations \omega [rad/s]')
% ylabel('PSD')
% figure;
% semilogy(w,abs(S_tt),'LineWidth',1.5); % Graphe PSD_thetatheta
% title( 'Nodal PSD S_{\theta\theta}(\omega)' )
% xlabel('Pulsations \omega [rad/s]')
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
    
    %D) Calcul de la vitesse du vent critique au flottement:
    if MAX_sigma_tt <= sigma_tt(j)
        MAX_sigma_tt = sigma_tt(j);
        U_crit = U(j); % 76.86 lorsque nous sommes précis
    end
    
    j = j + 1;
end

figure;
hold on;
plot(U,abs(sigma_zz),'-o','LineWidth',1.5);  % Graphe de l'écart-type sigma_zz
plot(U,abs(sigma_tt*B/2),'-o','LineWidth',1.5); % Graphe de l'écart-type sigma_thetatheta
xline(U_crit,'r','LineWidth',1.5);
title('Ecart-type de la réponse');
xlabel('Vitesse moyenne du vent U [m/s]');
ylabel('Ecart-type');
legend('\sigma_{zz}','\sigma_{\theta\theta} * B/2','U_{critique}')
hold off;
