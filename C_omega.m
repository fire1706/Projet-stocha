function C = C_omega(w,U,B)

a0 = 1;
a1 = -0.165;
a2 = -0.335;
b1 = 0.0455;
b2 = 0.3; 

w_1 = 2*U*b1/B;
w_2 = 2*U*b2/B;

C = a0 + a1*(1i*w)/(1i*w + w_1) + a2*(1i*w)/(1i*w + w_2); % Approximation de la fonction circulatoire de Theodorsen par Jones

end