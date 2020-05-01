function [PSD_zz,F_zz,PSD_tt,F_tt,PSD_ww,F_ww] = PSD_omega_2(time,omega_2,Zout,w)

dt = time(2)-time(1);

[PSD_z,F_z] = pwelch(Zout(:,1),[],[],length(omega_2),1/dt); % PSD de z(t)
F_zz = 2*pi*F_z; % Transformations de temporel vers fréquentiel pour la comparaison dans le domaine fréquentiel
PSD_zz = PSD_z/(4*pi);

[PSD_t,F_t] = pwelch(Zout(:,2),[],[],length(omega_2),1/dt); % PSD de theta(t)
F_tt = 2*pi*F_t; % Transformations de temporel vers fréquentiel pour la comparaison dans le domaine fréquentiel
PSD_tt = PSD_t/(4*pi);

[PSD_w,F_w] = pwelch(w,[],[],length(omega_2),1/dt);
F_ww = 2*pi*F_w;
PSD_ww = PSD_w/(4*pi);

end