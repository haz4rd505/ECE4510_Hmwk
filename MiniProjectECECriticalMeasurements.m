%State Estimation for the 5 bus system 
syms Vone Vtwo Vthree Vfour Vfive Dtwo Dthree Dfour Dfive

YbusG = [3 -1 0 0 -2; -1 3 0 -2 0; 0 0 0 0 0; 0 -2 0 3 -1; -2 0 0 -1 3];
YbusB = [-6 2 0 0 4; 2 -12.05 5.5 4 0; 0 5.5 -5 0 0; 0 4 0 -5.5 2; 4 0 0 2 -6];
Pvals = [0, -0.7, -0.4, 1, -0.6];
Qvals = [0, -0.2, -0.1, 0, -0.3];

%From Homework 5 we have the "measurements" we're going to be using for our
%Pinj,Qinj,V,Angle,and Lineflow from Powerworld

v1=1.05;
v2=.917;
v3=.988;
v4=1.05;
v5=.966;
a1=0;
a2=-7.43;
a3=-11.973;
a4=-0.059;
a5=-3.390;
Pinj1 = .817;Qinj1 = .411;Pinj2 = -.700;Qinj2 = -.200;Pinj3=-.400;Qinj3=-.100;
Pinj4 = 1.000;Qinj4 = -.146;Pinj5=-.600;Qinj5=-.300;
P12=.397;P21=.363;Q12=.171;Q21=-.103;
P24=.721;P42=.788;Q24=-.208;Q42=.341;
P23=.395;P32=.395;Q23=.122;Q32=-0.089;
P45=.208;P54=.197;Q45=.120;Q54=-.100;
P15=.420;P51=.399;Q15=.240;Q51=-0.198;

testvector = [a2,a3,a4,a5,v1,v2,v3,v4,v5];
testvector(1:4) = testvector(1:4)*(2*pi)/360


D = [0,  Dtwo, Dthree, Dfour, Dfive];
V = [Vone, Vtwo, Vthree, Vfour, Vfive];
% Setting up the equations for State Estimation.
%% Pinjected Powerflow equations
Pvec={};
F=0;
for kk = 1:5
    for ll = 1:5
        P = abs(V(kk))*abs(V(ll))*(YbusG(kk,ll)*cos(D(kk) - D(ll)) + YbusB(kk,ll)*sin(D(kk)-D(ll)));
        F = F + P ;
    end
    Pvec{kk} = F;
    F = 0;
end
%% Qinjected Powerflow Equations
Qvec = {};
G=0;
for ii = 1:5
    for jj = 1:5
        Q = abs(V(ii))*abs(V(jj))*(YbusG(ii,jj)*sin(D(ii) - D(jj)) - YbusB(ii,jj)*cos(D(ii)-D(jj)));
        G = G + Q ;
    end
    Qvec{ii} = G;
    G = 0; 
end

%%

%% Small g and b for only admittances, not the Ybus Matrix
smallG = YbusG;
for ii =1:5 
    smallG(ii,ii) = -smallG(ii,ii);
end
smallG = -smallG;

smallB = YbusB;
for ii =1:5 
    smallB(ii,ii) = -smallB(ii,ii);
end
smallB = -smallB;
%%
Pmatrix = zeros(5,5);
for kk = 1:5
    for ll = 1:5
    P = V(kk)^2*(smallG(kk,ll)) - V(kk)*V(ll)*(smallG(kk,ll)*cos(D(ll)-D(kk)) - smallB(kk,ll)*sin(D(ll)-D(kk)));
    Pmat(kk,ll) = P; 
    %Pmatrix(kk,ll)=subs(Pmat(kk,ll),[Dtwo,Dthree,Dfour,Dfive,Vone,Vtwo,Vthree,Vfour,Vfive],testvector)
    end
end
%%
Qmatrix = zeros(5,5);
for kk = 1:5
    for ll = 1:5
        Q = -V(kk)^2*(smallB(kk,ll)) + V(kk)*V(ll)*(-smallG(kk,ll)*sin(D(kk) - D(ll)) + smallB(kk,ll)*cos(D(kk)-D(ll)));
        Qmat(kk,ll) = Q;
       %Qmatrix(kk,ll)=subs(Qmat(kk,ll),[Dtwo,Dthree,Dfour,Dfive,Vone,Vtwo,Vthree,Vfour,Vfive],testvector)
    end 
end


%% With Q4,Q23 and Q32 removed
x = [0 0 0 0 1 1 1 1 1];
%Zwith error and without
errterm = normrnd(0,.005,1,16)
%Adding a single gross error
errterm(8) = 4;
%Z=[Pinj1;Pinj2;Pinj3;Pinj4;Pinj5;Qinj1;Qinj2;Qinj3;Qinj5;P12;-P21;;P23;-P32;-P24;P42;P45;-P54;P15;-P51;Q12;Q21;Q24;Q42;Q45;Q54;Q15;Q51;v1;v2;v3;v4;v5;]
Z=[Pinj1;Pinj4;Pinj5;Qinj1;Qinj5;P12;-P21;P23;-P24;P45;Q24;v1;v2;v3;v4;v5;]+transpose(errterm)
%+transpose(errterm)
%Q45;Q54

littleh = [Pvec{1},Pvec{4},Pvec{5},Qvec{1},Qvec{5},...
           Pmat(1,2),Pmat(2,1),Pmat(2,3),Pmat(2,4),Pmat(4,5),...
           Qmat(2,4),Vone,Vtwo,Vthree,Vfour,Vfive]
%Adding Error to the little hx
%Qmat(4,5),Qmat(5,4)

       
 
%hx = subs(littleh,[Dtwo,Dthree,Dfour,Dfive,Vone,Vtwo,Vthree,Vfour,Vfive],x)
bigH = jacobian(littleh,[Dtwo;Dthree;Dfour;Dfive;Vone;Vtwo;Vthree;Vfour;Vfive]);
%Hx = subs(bigH,[Dtwo,Dthree,Dfour,Dfive,Vone,Vtwo,Vthree,Vfour,Vfive],x)

Rmat = eye(length(Z))*(.005)^2;
Wmat = inv(Rmat)


%%
error = 1
x = [0 0 0 0 1 1 1 1 1];
i=0
while error > 1e-6
    hx = double(subs(littleh,[Dtwo,Dthree,Dfour,Dfive,Vone,Vtwo,Vthree,Vfour,Vfive],x))

    Hx = double(subs(bigH,[Dtwo,Dthree,Dfour,Dfive,Vone,Vtwo,Vthree,Vfour,Vfive],x));
    
    flip = double(inv(transpose(Hx)*Wmat*Hx))
    deltax = flip*double(transpose(Hx)*Wmat*(Z-transpose(hx)))
    x = x + transpose(deltax)
    error = norm(deltax)
    i = i +1
end
%%
Angles = transpose(x(1:4))*360/(2*pi);
Voltages = transpose(x(5:end));
StateEstimationAngles = table(Angles)
StateEstimationVoltages=table(Voltages)
[Pinj1;Pinj4;Pinj5;Qinj1;Qinj5;P12;-P21;P23;-P24;P45;Q24;v1;v2;v3;v4;v5;]
residuals = Z-transpose(hx)
figure(1)
bar(residuals)
set(gca,'XTick',[1:length(Z)],'XTickLabel',{'P1','P4','P5','Q1','Q5','P12','P21','P23','P24','P45','Q24','V1','V2','V3','V4','V5'})
title('Residuals of Critical Measurement with one Gross Error')

%%

%%Displaying the results
Zlabels = {'Pinj1','Pinj4','Pinj5','Qinj1','Qinj5','P12','-P21','P23','-P24','P45','Q24','v1','v2','v3','v4','v5'}
hxfinal = transpose(hx);
table(transpose(Zlabels),Z,hxfinal)
