function yita0=yita_dH_dG_mono(p,T,SA_mono,DMA,CS,dH,dG)
% cal pre-param
kb=1.381e-23;
beta=1.126e-15*(T/298.15).^0.5;
gamma=cal_evap_rate_dH_dG(p,298.15,dH,dG)*(T/298.15).^(0.5).*exp(dH*4185.85/8.314./T-dH*4185.85/8.314/298.15);

% unit normalization
DMA=DMA*1e-12*p/kb./T;
CS=1.3*CS;

% cal med-param
rCS=CS./beta;
rGA=gamma./beta;
% cal J
A=0.96*DMA+rGA+0.63*rCS;
B=0.86*SA_mono-rGA-0.63*rCS;
C=-0.86*SA_mono;
T=sqrt(B.^2-4*A.*C);
yita0=1-(T-B)./A/2;
end