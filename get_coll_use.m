% Collision coefficients
function K=get_coll_use(p,T)
K = nan(24,24);

ai=[1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4];
bi=[0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4];
for i=1:24
    for j=1:24
        K(i,j)=calbeta(p,T,ai(i),bi(i),ai(j),bi(j));
    end
end

%K=K*(T/T0)^(1/2);
end
