function Diff=calDiff(p,T,dp)
kB=1.381e-23;
Cs=calCs(p,T,dp);
vis=calvis(T);
Diff=kB*T.*Cs./(3*pi*vis*dp);



end