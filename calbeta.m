function beta=calbeta(p,T,a1,b1,a2,b2)

[m1,dp1]=calm_dp(a1,b1);
[m2,dp2]=calm_dp(a2,b2);
Diff1=calDiff(p,T,dp1);
Diff2=calDiff(p,T,dp2);
v1=calvelocity(m1,T);
v2=calvelocity(m2,T);
l1 = 8*Diff1/pi./v1;
l2 = 8*Diff1/pi./v2;
g1 = sqrt(2)./(3*dp1*l1).* ( (dp1+l1).^3 - (dp1^2+l1.^2).^(3/2) ) -dp1;
g2 = sqrt(2)./(3*dp2*l2).* ( (dp2+l2).^3 - (dp2^2+l2.^2).^(3/2) ) -dp2;
a = 2*pi*(Diff1+Diff2)*(dp1+dp2);
b = (dp1+dp2)./(dp1+dp2+2*sqrt(g1.^2+g2.^2));
c = 8*(Diff1+Diff2)./sqrt(v1.^2+v2.^2)/(dp1+dp2);
beta=2.3*a./(b+c);

end