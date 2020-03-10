function S = S_omega(w)
U = 0:1:100;
L_w = 20;
I_w = 0.05;
sigma = U*I_w;

for i=1:1:101
    jean(i) = (w*L_w)/(2*pi*U(i));
    A(i) = sigma(i)*sigma(i)*L_w/pi/U(i);
    B(i) = 1+ 755.2*jean*jean;
    C(i) = 1+283.2*jean*jean;
    S(i) = A*B/C^(11/6);
end
end

