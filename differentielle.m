function [dz] = differentielle(t,z,A_final,M_final,f_b)
dz = zeros(10,1);
b0 = 0;
b1 = 0.0455;
b2 = 0.3;
b = [b0,b1,b2];
t_max = 600;
N = t_max;
time = linspace(0,t_max,N);
dt = time(2)-time(1);
%f_b_2 = interp1(time,f_b,t);
i     = floor(t/dt); % ici on repère la position de t dans le vecteur time
fb  = interp1(time(1+[i i+1]),f_b(:,1+[i i+1]).',t).';
%f_b_2  = transpose(interp1(time(1+[i i+1]),permute(f_b(:,:,1+[i i+1]),[2,1,3]),t)); 
dz(:,1) = (A_final(:,1)*z(1) + fb(:,1))/(M_final(:,1));

end
