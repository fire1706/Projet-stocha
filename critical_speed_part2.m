function [U_crit,sigma_z,sigma_theta] = critical_speed_part2(VALEURS_PROPRES,U,sigma_z,sigma_theta)

% U = 10:1:100;
% %U_crit = 76.86; % Valeur calcul�e dans la partie 1
% sigma_z(68:length(U)) = 0; % Etant donn� que le syst�me devient instable pass� la vitesse critique moyenne du vent U_crit, on va d�finir � z�ro les valeurs d'�cart-type pour U > U_crit
% sigma_theta(68:length(U)) = 0;
k = 2;
while k <= length(U)
    if VALEURS_PROPRES(1,k-1)-VALEURS_PROPRES(1,k) <  VALEURS_PROPRES(1,k-1)
        critical_position = k;
        U_crit = U(critical_position);
        break;
    end
    k = k + 1;
end
sigma_z(critical_position:length(U)) = 0; % Etant donn� que le syst�me devient instable pass� la vitesse critique moyenne du vent U_crit, on va d�finir � z�ro les valeurs d'�cart-type pour U > U_crit
sigma_theta(critical_position:length(U)) = 0;
end
