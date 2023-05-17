function [J_lyy,J1,J2,dt1,dt2,dt3,Conc1,Conc2]= get_J_comp(cm3A,pptB,CS,temp_conv,dH_lyy,dG_lyy,dG_tot,num1,num2)
if dH_lyy~=-24.82||dG_lyy~=-13.54
    error('dH or dG need to be updated!')
end
kB=1.381e-23;


cd('D:\researchwork\I_into_WRFChem');
Beta=get_coll_use(101325,temp_conv);

[Rel,dG_rel]=get_rel_dG(dG_tot);
Gamma=get_evap_use(dG_rel,Beta,temp_conv);
Conc1=zeros(1,24);
Conc1(1)=cm3A*1e6;
Conc1(5)=pptB*1e-12*101325/kB./temp_conv;
Conc2=zeros(1,7);
Conc2(1)=cm3A*1e6;
Conc2(2)=pptB*1e-12*101325/kB./temp_conv;
Gamma2=zeros(7,7);
Gamma2(1,2)=cal_evap_rate_dH_dG(101325, temp_conv,dH_lyy,dG_lyy);
Gamma2(2,1)=cal_evap_rate_dH_dG(101325, temp_conv,dH_lyy,dG_lyy);
Beta2=get_coll_use_dyn(101325,temp_conv);
Rel2=nan(7,7);
Rel2(1,2)=3;
Rel2(2,1)=3;
Rel2(1,3)=4;
Rel2(3,1)=4;
Rel2(4,2)=5;
Rel2(2,4)=5;
Rel2(3,5)=6;
Rel2(5,3)=6;
Rel2(5,5)=7;
Rel2(6,3)=7;
Rel2(3,6)=7;
[Ai,Di]=calAD(1:length(Beta));
[~,Dp_tot]=Bulk_m_dp(Ai,Di);
ai=[1,0,1,2,2,3,4];
bi=[0,1,1,1,2,3,4];
[~,Dp_tot2]=Bulk_m_dp(ai,bi);
TI=1;
Conc1_rec=Conc1;
Conc2_rec=Conc2;
t1=cputime;

for i=1:num1
Conc_new=Manual_Dynamics(TI,Conc1_rec,Gamma,CS,Dp_tot,Rel,Beta);
Conc_new(1)=nansum(Conc1_rec(1:5:21))-nansum(Conc_new(6:5:21));
Conc_new(5)=Conc1_rec(5);
Conc1_rec=Conc_new;
end
t2=cputime;
dt1=t2-t1;
t1=cputime;
for i=1:num2
Conc_new2=Manual_Dynamics2(TI,Conc2_rec,Gamma2,CS,Dp_tot2,Rel2,Beta2);
Conc_new2(1)=nansum(Conc2_rec([1,3]))-nansum(Conc_new2(3));
Conc_new2(2)=Conc2_rec(2);
Conc2_rec=Conc_new2;
end
t2=cputime;
dt2=t2-t1;
Conc1=Conc1_rec;
Conc2=Conc2_rec;
%%
Coll1=Conc1'*Conc1.*Beta;
Coll2=Conc2'*Conc2.*Beta2;
J1=Coll1(12,12)/2-Gamma(12,12)*Conc1(24);
for j=1:11
    J1=J1+Coll1(j,24-j)-Gamma(j,24-j)*Conc1(24);
end
J2=Coll2(5,5)/2+Coll2(3,6);

%%
t1=cputime;
J_lyy=JLYYnew_dH_dG(101325,temp_conv,cm3A*1e6,pptB,CS,dH_lyy,dG_lyy);
t2=cputime;
dt3=t2-t1;
end