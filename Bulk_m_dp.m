function [M_tot,Dp_tot]=Bulk_m_dp(Ai,Di)
    M_tot=nan(length(Ai),1);
    for i=1:length(Ai)
        [M_tot(i),Dp_tot(i)]=calm_dp(Ai(i),Di(i));
    end
end