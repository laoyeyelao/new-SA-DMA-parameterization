function [Rel,dG_rel]=get_rel_dG(dG)
    [Ai,Di]=calAD(1:length(dG));
    maxAD=sqrt(length(dG)+1)-1;
    dG_rel=nan(length(dG),length(dG));
    Rel=nan(length(dG),length(dG));
    for i=1:length(dG)
        for j=1:length(dG)
            nA=Ai(i)+Ai(j);
            nD=Di(i)+Di(j);
            if nA<=maxAD&&nD<=nA
                Rel(i,j)=nD*(maxAD+1)+nA;
                dG_rel(i,j)=dG(nD*(maxAD+1)+nA)-dG(i)-dG(j);
            end
        end
    end

end