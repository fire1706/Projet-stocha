function [H_zz,H_tt] = FRF_part2(omega_2,A_final,M_final)

H = inv(1i*omega_2*M_final - A_final); % FCT de transfert
H_zz = H(1,1);
H_tt = H(2,2);

end
