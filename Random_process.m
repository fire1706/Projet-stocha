function [w,S_w,W] = Random_process(d_omega,omega_2,N,U,L_w,I_w)

k = 1;
while k <= length(omega_2)
theta(k) = 2*pi*rand(1,1);
S_w(k) = S_omega(omega_2(k),U,L_w,I_w); % PSD de Von Karman en fréquentiel 
W(k) = N*sqrt(S_w(k)*d_omega)*exp(1i*theta(k)); % W(omega) = composante fluctuante du vent en fréquentiel
k = k + 1;
end

if (N/2 - round(N/2)) == 0 
    W(1) = 0; % Conditions (pour vérifier la symétrie) pour la décomposition en série de Fourier d’une fonction continue: N est pair
    W(N/2 + 1) = 0;
    W(N/2 + 2 : N)  = conj(W(N/2:-1:2));
else
    W(1) = 0; % Conditions (pour vérifier la symétrie) pour la décomposition en série de Fourier d’une fonction continue: N est impair
    W((N+3)/2 : N)  = conj(W((N+1)/2:-1:2)); 
end

w = ifft(W,'symmetric'); % w(t) = composante fluctuante du vent en temporel

end