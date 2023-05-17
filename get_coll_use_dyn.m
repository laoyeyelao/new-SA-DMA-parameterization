% Collision coefficients
function K=get_coll_use_dyn(p,T)
K = nan(7,7);

ai=[1,0,1,2,2,3,4];
bi=[0,1,1,1,2,3,4];
for i=1:7
    for j=1:7
        K(i,j)=calbeta(p,T,ai(i),bi(i),ai(j),bi(j));
    end
end

%K=K*(T/T0)^(1/2);
end
