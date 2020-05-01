function [t,Zout,VALEURS_PROPRES,A_final,M_final] = Diff_equation_for_ode45(omega_2,time,U,B,rho,W,w,M,K,C,Flow)

% F_b_omega = 1/4*rho*U*B*[4*pi*W.', pi*B*W.']; % Force de turbulence en fréquentiel
% f_b_omega2 = ifft(F_b_omega,'symmetric'); % Devrait être égal à f_b_omega car opérations linéaires => OK!
f_b_omega = 1/4*rho*U*B*[4*pi*w.', pi*B*w.'];% Force de turbulence en temporel
f_b_omega = transpose(f_b_omega); %On transpose f_b_omega pour pouvoir extraire L_b et M_b comme étant des vecteurs lignes
L_b = f_b_omega(1,:);
M_b = f_b_omega(2,:);

if Flow == 1
    a0 = 1;
    a1 = 0;
    a2 = 0;
    b0 = 0;
    b1 = 0;
    b2 = 0;
elseif Flow == 0
    a0 = 1;
    a1 = -0.165;
    a2 = -0.335;
    b0 = 0;
    b1 = 0.0455;
    b2 = 0.3;
end

Mat_B = -2*U/B *diag([b0 b1 b2 b0 b1 b2]);
Mat_A =[0,a0,0,0;0,a1,0,0;0,a2,0,0;0,0,-a0/U,a0*B/(4*U);0,0,-a1/U,a1*B/(4*U);0,0,-a2/U,a2*B/(4*U)];

q = pi*rho*U^2*B;
Mat_S = q*[0,0,0,0,0,0;0,0,0,0,0,0;1,1,1,1,1,1;B/4,B/4,B/4,B/4,B/4,B/4];

M_conj = [pi*rho*B^2/4,0;0,pi*rho*B^4/128];
C_conj = [0,-q*B/(4*U);0,q*B^2/(16*U)];

M_0 = [eye(2),zeros(2,2);zeros(2,2),M+M_conj];
A_0 = [zeros(2,2),eye(2);-K,-(C+C_conj)];

M_final = [M_0,zeros(4,6);-Mat_A,eye(6)];
A_final = [A_0,Mat_S;zeros(6,4),Mat_B];
VALEURS_PROPRES = real(eig(M_final\A_final)); % Pour la vérification que les parties réelles des valeurs propres soient négatives

f_b = zeros(10,length(omega_2));
f_b([3 4],:) = [L_b;M_b]; % Forces de turbulence en temporelle 

tspan = time;
z0 = zeros(10,1); % structure initialement immobile
[t,Zout] = ode45(@(t,Zout) differentielle(t,Zout,A_final,M_final,f_b,time) ,tspan,z0); %ODE45 résout notre équation différentielle " M*z_dot = A*z + f_b"

end