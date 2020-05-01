function [U_crit] = critical_speed(sigma_tt,sigma_zz,U)
U_crit = 0;

cumul_MAXtt = cummax(sigma_tt);
cumul_MAXzz = cummax(sigma_zz);
i = 2;
while i <= length(U)-1
    if cumul_MAXtt(i-1) == cumul_MAXtt(i) && cumul_MAXzz(i-1) == cumul_MAXzz(i)
        U_crit = U(i-1);
        break;
    end
    i = i + 1;
end

end