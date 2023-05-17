function lam=calmfp(p,T)
p0=101325;
T0=293.15;
lam=66.5e-9.*(p0./p).*(T/T0).^2.*((T0+110.4)./(T+110.4));



end