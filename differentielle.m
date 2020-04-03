function [dz] = differentielle(t,z,A_final,M_final,fb)
dz = zeros(10,1);
% t_max = 600;
% N = t_max;
% time = linspace(0,t_max,N);
% dt = time(2)-time(1);
% fb = interp1(time,f_b.',t).';
% i     = floor(t/dt); % ici on repère la position de t dans le vecteur time
%fb  = interp1(time(1+[i i+1]),f_b(:,1+[i i+1]).',t).';
%f_b_2  = transpose(interp1(time(1+[i i+1]),permute(f_b(:,:,1+[i i+1]),[2,1,3]),t)); 
dz(:,1) = M_final(:,:)\(A_final(:,:)*z(:) + fb(:,1));

end
