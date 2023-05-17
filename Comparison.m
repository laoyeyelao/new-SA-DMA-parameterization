function [J1_par,J1_ACDC,J1_dynamic]=Comparison(T,SA,DMA,CS)
dH=-24.82;
dG=-13.54;
J1_par=zeros(length(T),1);
J1_ACDC=zeros(length(T),1);
J1_dynamic=zeros(length(T),1);
dt1_ACDC=zeros(length(T),1);
dt2_dynamic=zeros(length(T),1);
Conc1=cell(length(T),1);
for k=1:length(T)
dG_tot=get_dG(T(k));
[J1_par(k),J1_ACDC(k),J1_dynamic(k),dt1_ACDC(k),dt2_dynamic(k),~,Conc1{k},~]= get_J_comp(SA(k),DMA(k),CS(k),T(k),dH,dG,dG_tot,500000,500000);
end
end