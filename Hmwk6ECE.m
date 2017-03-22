% Solving Power Flow Equations from Hmwk5
%%
YbusG = [3 -1 0 0 -2; -1 3 0 -2 0; 0 0 0 0 0; 0 -2 0 3 -1; -2 0 0 -1 3];
YbusB = [-6 2 0 0 4; 2 -12.05 5.5 4 0; 0 5.5 -5 0 0; 0 4 0 -5.5 2; 4 0 0 2 -6];
Pvals = [0, -0.7, -0.4, 1, -0.6];
Qvals = [0, -0.2, -0.1, 0, -0.3];
% Note, because we defined the variables this script can only be used for a
% 5 bus system with one slack bus, one pv bus, and 3 pq busses.
syms Vtwo Vthree Vfive Dtwo Dthree Dfour Dfive

% Writing in the powerflow equations
%%

D = [0, Dtwo, Dthree, Dfour, Dfive];
V = [1.05, Vtwo, Vthree, 1.05, Vfive];

%%
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


jcobi2 = jacobian([Pvec{2},Pvec{3},Pvec{4},Pvec{5},Qvec{2},Qvec{3},Qvec{5}],[Dtwo;Dthree;Dfour;Dfive;Vtwo;Vthree;Vfive]);
%initial angles
done = 0;
dtwo = 0;
dthree = 0;
dfour= 0;
dfive = 0;
%initial voltages
vone = 1.05
vtwo = 1
vthree = 1
vfour = 1.05
vfive = 1
initialguess = [dtwo;dthree;dfour;dfive;vtwo;vthree;vfive]
%% 
jcobivals = zeros(7,7)
err =10
i = 0

while err > 1e-16
    


for pp = 1:7 
    for qq = 1:7  
jcobivals(pp,qq) =  subs(jcobi2(pp,qq),[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess));
    end
end
invertedjcobi = inv(jcobivals)
fofx = [subs(Pvec{2},[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess));
        subs(Pvec{3},[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess));
        subs(Pvec{4},[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess));
        subs(Pvec{5},[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess));
        subs(Qvec{2},[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess));
        subs(Qvec{3},[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess));
        subs(Qvec{5},[Dtwo,Dthree,Dfour,Dfive,Vtwo,Vthree,Vfive],transpose(initialguess))];
    
    deltax = double(invertedjcobi*fofx);
    initialguess = initialguess - deltax
i = i+1
    err = abs(deltax(7))
    i = i+1
end

initialguess(1:4) = initialguess(1:4)*360/(2*pi);
%initialguess(1:4) = initialguess(1:4)*2*pi/360
finaloutput = initialguess



