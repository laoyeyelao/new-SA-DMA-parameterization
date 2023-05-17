function Evap=get_evap_use(dG_rel,Beta,T)
    Rg=8.314;
    NA=6.02e23;
    Evap=Beta.* (101325/Rg./298.15) * NA.*exp(4185.85*dG_rel/Rg/T);
end