%% Partie 1 => fréquentielle
clear all;
%DONNEES:
%1) Propriétés PONT
i = 1;
w = 0:0.01:2;
m = 22740;
J = 2.47*10^6;
B = 31;
xhi = 0.003;
f_z = 0.10;
f_theta = 0.278;

%2) Propriétés VENT
rho = 1.22;
L_w = 20;
U = 0:1:100;%U = [0,100];
j = 1;
I_w = 0.05;

%A) Fct de transfert:
w_0 = [2*pi*f_z,2*pi*f_z;2*pi*f_theta,2*pi*f_theta];
M = [m,0;0,J];
K = w_0.^2.*M;
C = 2*xhi*sqrt(K*M);%C = 2*xhi*w_0.*M;

while j <= length(U)
    
    % M*x_ddot + C*x_dot + K*x = f
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
        C_w(i) = C_omega(w(i),U(j),B);
        q = pi*rho*U(j)^2*B;
        Grosse_matrice_desesmorts = [B*w(i)^2/(4*U(j)^2) - C_w(i)*1i*w(i)/U(j) , B/(4*U(j))*1i*w(i) + C_w(i)*(1+B/(4*U(j))*1i*w(i)) ; -B/4*C_w(i)*1i*w(i)/U(j) , B/4*(B^2*w(i)^2/(32*U(j)^2) - B/(4*U(j))*1i*w(i) +  C_w(i)*(1+B/(4*U(j))*1i*w(i)))];
        H(:,:,i) = inv(-w(i)^2*M + 1i*w(i)*C + K - q*Grosse_matrice_desesmorts);
        H_z1(i) = H(1,1,i);
        H_t2(i) = H(2,2,i);
%       H_t1(i) = H(2,1,i);
%       H_z2(i) = H(1,2,i);
        %F_se_omega = q*Grosse_matrice_desesmorts*[POINT_interro1;POINT_interro2];
        % F_b_omega = 1/4*rho*U*B*[4*pi;pi*B]*W_omega;
        % H = X/fb;
        
        %B) Calcul densité de puissance
        S_w(i) = S_omega(w(i),U(j));
        F_b_omega = (1/4)*rho*U(j)*B*[4*pi;pi*B];
        S_b(:,:,i) = F_b_omega*S_w(i)*conj(transpose(F_b_omega));
        S_x(:,:,i) = H(:,:,i)*S_b(:,:,i)*conj(transpose(H(:,:,i)));
        S_zz(i) = S_x(1,1,i);
        S_tt(i) = S_x(2,2,i);
        
%         S_tz(i) = S_x(2,1,i);
%         S_zt(i) = S_x(1,2,i);
        i = i + 1;
    end
    % H_zz = sqrt(H_z1.^2+H_z2.^2);
    % H_tt = sqrt(H_t1.^2+H_t2.^2);
    
%     figure;
%     plot(w,abs(H_z1),'LineWidth',1.5);
%     title( 'Nodal FRF H_{zz}(\omega)' )
%     xlabel('Pulsations \omega [rad/s]')
%     ylabel('FRF')
%     figure;
%     plot(w,abs(H_t2),'LineWidth',1.5);
%     title( 'Nodal FRF H_{\theta\theta}(\omega)' )
%     xlabel('Pulsations \omega [rad/s]')
%     ylabel('FRF')
%     figure;
%     semilogy(w,abs(S_zz),'LineWidth',1.5);
%     title( 'Nodal PSD S_{zz}(\omega)' )
%     xlabel('Pulsations \omega [rad/s]')
%     ylabel('PSD')
%     figure;
%     semilogy(w,abs(S_tt),'LineWidth',1.5);
%     title( 'Nodal PSD S_{\theta\theta}(\omega)' )
%     xlabel('Pulsations \omega [rad/s]')
%     ylabel('PSD')
    sigma_zz(j) = sqrt(trapz(w,S_zz));
    sigma_tt(j) = sqrt(trapz(w,S_tt));
    
    %C) Vérification Theodorsen
%     figure;
%     hold on
%     plot(w,imag(C_w),'LineWidth',1.5);
%     plot(w,real(C_w),'LineWidth',1.5);
%     hold off
    
    j = j + 1;
    
end
grid;
% figure;
hold on;
plot(U,abs(sigma_zz),'-o','LineWidth',1.5);
plot(U,abs(sigma_tt*B/2),'-o','LineWidth',1.5);
xline(77,'r','LineWidth',1.5);
title( 'Standard deviation of response' );
xlabel('U [m/s]');
ylabel('Deviation');
hold off;
