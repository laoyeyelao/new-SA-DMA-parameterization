function J_lyy=JLYYnew_dH_dG(p,T,SA,DMA,CS,dH,dG)
% cal pre-param
kb=1.381e-23;
beta=1.126e-15*(T/298.15).^0.5;
gamma=cal_evap_rate_dH_dG(p,298.15,dH,dG)*(T/298.15).^(0.5).*exp(dH*4185.85/8.314./T-dH*4185.85/8.314/298.15);

% unit normalization
DMA=DMA*1e-12*p/kb./T;
CS=1.3*CS;

% cal med-param
rCS=CS./beta./SA;
rGA=gamma./beta./SA;
rDMA=DMA./(DMA+0.39*rCS.*SA);
% cal J
yita0=0.96*DMA./SA./(0.96*DMA./SA+rGA+0.86+0.63*rCS);
% A=0.45-0.77*rDMA+(0.38*rCS).^2;
% B=0.38*rCS+0.77*rDMA;
% pA=1.11*(((1+yita0./(yita0+0.28*rCS)).*((2*A.*yita0+B)/2./sqrt(A.*yita0.^2+B.*yita0+0.25)-0.38*rCS)-0.28*rCS./(yita0+0.28*rCS).^2.*(sqrt(A.*yita0.^2+B.*yita0+0.25)-0.38*rCS.*yita0-0.25)))+0.14;
% pB=-pA.*yita0+1.11*(2*yita0+0.28*rCS)./(yita0+0.28*rCS).*(sqrt(A.*yita0.^2+B.*yita0+0.25)-0.38*rCS.*yita0-0.5)+0.14*yita0+0.96*DMA./SA+rGA+0.86+0.63*rCS;
% yita=(-pB+(pB.^2+4*0.96*pA.*DMA./SA).^0.5)/2./pA;
AB=yita0.*SA;
theta=1+2*DMA./(1.16*DMA+0.46*rCS.*SA).*(SA-AB)./AB;
theta2=theta.*(2.22*AB+0.86*rCS.*SA)./(((1.11*AB+0.43*rCS.*SA).^2+1.12*theta.*AB.^2).^0.5+1.11*AB+0.43*rCS.*SA);
J_lyy=beta.*theta2.*AB.^4/2./(AB+0.39*rCS.*SA).*(0.23*theta2./(AB+0.39*rCS.*SA)+1./(AB+0.31*rCS.*SA));
end