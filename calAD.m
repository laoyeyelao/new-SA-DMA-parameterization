function [Ai,Di]=calAD(I)
[a,b]=size(I);
if min(a,b)~=1
    error('not proper I type!');
end
Ai=nan(length(I),1);
Di=nan(length(I),1);
ADnum=sqrt(length(I)+1);
for i=1:length(I)
    Ai(i)=rem(I(i),ADnum);
    Di(i)=floor(I(i)/ADnum);
end
end