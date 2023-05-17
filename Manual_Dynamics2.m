function Conc_new=Manual_Dynamics2(TI,Conc,Gamma,CS,Dp_tot,Rel,Beta)
    if TI~=0
    Coll=Conc'*Conc.*Beta;
    Evap_tot=zeros(length(Conc),1);
    CoagS=1.3*CS*(Dp_tot/Dp_tot(1)).^(-1.7);
    CoagS_tot=-CoagS';
    for i=1:length(Conc)
        for j=1:i
            if ~isnan(Rel(i,j))
                Evap_tot(Rel(i,j))=Evap_tot(Rel(i,j))-Gamma(i,j);
            end
        end
    end
    Conc_new=Conc.*exp((CoagS_tot+Evap_tot)*TI/2)';
    F_evap=(1-exp((CoagS_tot+Evap_tot)*TI/2))./(CoagS_tot+Evap_tot)/TI*2;
    for i=1:length(Conc_new)
        for j=1:i
            if ~isnan(Rel(i,j))
                    Conc_new(j)=Conc_new(j)+Gamma(i,j)*Conc_new(Rel(i,j))*F_evap(Rel(i,j))*TI;
                    Conc_new(i)=Conc_new(i)+Gamma(i,j)*Conc_new(Rel(i,j))*F_evap(Rel(i,j))*TI;
            end
        end
    end
    for i=1:length(Conc_new)
        for j=1:i
            if ~isnan(Rel(i,j))
                if i~=j
                    Conc_new(i)=Conc_new(i)-Coll(i,j)*TI;
                    Conc_new(j)=Conc_new(j)-Coll(i,j)*TI;
                    Conc_new(Rel(i,j))=Conc_new(Rel(i,j))+Coll(i,j)*TI;
                else
                    Conc_new(j)=Conc_new(j)-Coll(i,j)*TI;
                    Conc_new(Rel(i,j))=Conc_new(Rel(i,j))+Coll(i,j)/2*TI;
                end
            end
        end
    end
    Conc_new=Conc_new.*exp((CoagS_tot+Evap_tot)*TI/2)';
    Conc_new(Conc_new<0)=0;
    else
        Conc_new=Conc;
    end
end