function [dz] = differentielle(t,z,A_final,M_final,f_b,time)
dz = zeros(10,1);


%fb_t = interp1(time,f_b.',t).'; % On calcule par interpolation les valeurs de f_b à la position de t (donné par ODE45) au lieu du vecteur time

dt = time(2)-time(1);
i     = floor(t/dt); % ici on repère la position de t dans le vecteur time
if t == time(end)
    fb_t = f_b(:,end);
else
    fb_t = interp1(time(1+[i i+1]),f_b(:,1+[i i+1]).',t).';
end

dz(:,1) = M_final(:,:)\(A_final(:,:)*z(:) + fb_t(:,1)); % Equation différentielle du mouvement

end
