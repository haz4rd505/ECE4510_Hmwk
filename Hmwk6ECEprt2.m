
%% Calculting Pline for all lines
dtwo = finaloutput(1);
dthree = finaloutput(2);
dfour = finaloutput(3);
dfive = finaloutput(4);
vtwo = finaloutput(5);
vthree = finaloutput(6);
vfive = finaloutput(7);
initialguess(1:4) = initialguess(1:4)*2*pi/360
inradians = initialguess



Pinjone = subs(Pvec{1},[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(inradians));
Qinjone = subs(Qvec{1},[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(inradians));
Qinjfour = subs(Qvec{4},[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(inradians));

shunt = 1.05^2/2
Pinj = [double(Pinjone), Pvals(2:5)];
Qinj = [double(Qinjone), Qvals(2:3), double(Qinjfour), Qvals(5), shunt];

Real_losses = sum(Pinj)
Reactive_losses = sum(Qinj)

%%



Pmatrix = zeros(5,5);
for kk = 1:5
    for ll = 1:5
    P = V(kk)*V(ll)*(YbusG(kk,ll)*cos(D(kk) - D(ll)) + YbusB(kk,ll)*sin(D(kk)-D(ll)));
    Pmatrix(kk,ll)=subs(P,[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(inradians))
    end
end
%%
Ploss12 = Pmatrix(1,2) + Pmatrix(2,1)
Ploss23 = Pmatrix(2,3) + Pmatrix(3,2)
Ploss24 = Pmatrix(2,4) + Pmatrix(4,2)
Ploss15 = Pmatrix(1,5) + Pmatrix(5,1)
Ploss54 = Pmatrix(5,4) + Pmatrix(4,5)
sum([Ploss12,Ploss23,Ploss24,Ploss15,Ploss54])

Qmatrix = zeros(5,5);
for ii = 1:5
    for jj = 1:5
        Q = V(ii)*V(jj)*(YbusG(ii,jj)*sin(D(ii) - D(jj)) - YbusB(ii,jj)*cos(D(ii)-D(jj)));
        Qmatrix(ii,jj)=subs(Q,[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(inradians));
    end 
end
%%
Ssquared = Qmatrix.^2 + Pmatrix.^2

Vvector = [vone, vtwo,vthree,vfour,vfive]
I12=(vone-vtwo)*(YbusG(1,2)^2+YbusB(1,2)^2)^(1/2)
I15=(vone-vfive)*(YbusG(1,5)^2+YbusB(1,5)^2)^(1/2)
I24=(vtwo-vfour)*(YbusG(2,4)^2+YbusB(2,4)^2)^(1/2)
I45=(vfour-vfive)*(YbusG(4,5)^2+YbusB(4,5)^2)^(1/2)
I23=(vthree-vtwo/.9091)*YbusB(2,3)