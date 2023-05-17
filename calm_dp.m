function [m,dp]=calm_dp(a,b)
NA=6.02e23;
rhoi=[1830,680]; % kg/m3
mi=[0.098,0.045]/NA; %kg/mol
num=[a,b];
vi=mi./rhoi;% m3
m=sum(mi.*num);
v=sum(vi.*num);
dp=(6*v/pi)^(1/3);




end