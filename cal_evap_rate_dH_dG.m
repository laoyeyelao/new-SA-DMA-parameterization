function gamma=cal_evap_rate_dH_dG(p, T,dH,dG)
    dfGm11=4185.85*(dG*(T/298.15)+(1-T/298.15)*(dH));
    Rg=8.314;
    NA=6.02e23;
    beta=calbeta(p,T,1,0,0,1);
    gamma = beta .* (101325/Rg./298.15) * NA.* exp(dfGm11/Rg./T); % add NA here because beta and gamma is defined based on molecule number
end