%%Fast Decoupled Load Flow
YbusG = [3 -1 0 0 -2; -1 3 0 -2 0; 0 0 0 0 0; 0 -2 0 3 -1; -2 0 0 -1 3];
YbusB = [-6 2 0 0 4; 2 -12.05 5.5 4 0; 0 5.5 -5 0 0; 0 4 0 -5.5 2; 4 0 0 2 -6];
Bprime = YbusB(2:5,2:5);
invertedBprime = inv(Bprime);
Bdprime = [-12.05 5.5 0; 5.5 -5 0; 0 0 -6];
invertedBdprime = inv(Bdprime);
Pvals = [0, -0.7, -0.4, 1, -0.6];
Qvals = [0, -0.2, -0.1, 0, -0.3];

syms Vtwo Vthree Vfive Dtwo Dthree Dfour Dfive

D = [0, Dtwo, Dthree, Dfour, Dfive];
V = [1.05, Vtwo, Vthree, 1.05, Vfive];

G = 0;
Qvec = {};
for ii = 1:5
    for jj = 1:5
        Q = V(ii)*V(jj)*(YbusG(ii,jj)*sin(D(ii) - D(jj)) - YbusB(ii,jj)*cos(D(ii)-D(jj)));
        G = G + Q ;
    end
    Qvec{ii} = G - Qvals(ii);
    G = 0; 
end

%%
F=0;
Pvec={};
for kk = 1:5
    for ll = 1:5
        P = V(kk)*V(ll)*(YbusG(kk,ll)*cos(D(kk) - D(ll)) + YbusB(kk,ll)*sin(D(kk)-D(ll)));
        F = F + P ;
    end
    Pvec{kk} = F - Pvals(kk);
    F = 0;
end
%%
%initial angles
done = 0;
dtwo = 0;
dthree = 0;
dfour= 0;
dfive = 0;
%initial voltages
vone = 1.05;
vtwo = 1;
vthree = 1;
vfour = 1.05;
vfive = 1;

initialguess = [dtwo;dthree;dfour;dfive;vtwo;vthree;vfive];
err = 10;
i=1;
%%
while err >1e-6
    deltaPoV = [subs((Pvec{2}/Vtwo),[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess));
        subs((Pvec{3}/Vthree),[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess));
        subs((Pvec{4}/vfour),[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess));
        subs((Pvec{5}/Vfive),[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess))];
    deltaQoV =[ subs((Qvec{2}/Vtwo),[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess));
        subs((Qvec{3}/Vthree),[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess));
        subs((Qvec{5}/Vfive),[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess))];
    
    deltaangle=double(invertedBprime*deltaPoV);
    deltavoltage=double(invertedBdprime*deltaQoV);
    
    delta = [deltaangle;deltavoltage];
    initialguess = double(initialguess + delta)
    err = abs(norm(delta));
    
end

initialguess(1:4) = initialguess(1:4)*360/(2*pi);
finaloutput = initialguess
%%



















