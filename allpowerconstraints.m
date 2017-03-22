function [c,ceq] = allpowerconstraints(x)

 y11 = 2*.00744/(.00744^2+.0372^2) - 2*.0372i/(.00744^2+.0372^2)+.0775i;
 y12 = -.00744/(.00744^2+.0372^2) + .0372i/(.00744^2+.0372^2);
 y21 = y12;
 y31 = y21;
 y13 = y31;
 y22 = (.00744/(.00744^2+.0372^2)+.01272/(.0636^2+.01272^2)) - (.0636/(.0636^2+.01272^2) + .0372/(.00744^2+.0372^2) - (.1275+.0775)/2)*i;
 y23 = -.01272/(.0636^2+.01272^2) + .0636i/(.0636^2+.01272^2);
 y32 = y23;
 y33 = -y32 - y31 + (.1275+.0775)/2*i;
Ybus = [ y11,y12,y13;y21,y22,y23;y31,y32,y33];
YbusB = imag(Ybus);
YbusG = real(Ybus);

V = [x(5), x(6),x(7)];
D = [0,x(8),x(9)];

Pvec={};
F=0;

%% Real Power Injected at a Bus
for ii = 1:3
    for jj = 1:3
        P = V(ii)*V(jj)*(YbusG(ii,jj)*cos(D(ii) - D(jj)) + YbusB(ii,jj)*sin(D(ii)-D(jj)));
        F = F + P ;
    end
    Pvec{ii} = F;
    F = 0;
end
%% Reactive Power Injected at a Bus
Qvec = {};
G=0;
for ii = 1:3
    for jj = 1:3
        Q = V(ii)*V(jj)*(YbusG(ii,jj)*sin(D(ii) - D(jj)) - YbusB(ii,jj)*cos(D(ii)-D(jj)));
        G = G + Q ;
    end
    Qvec{ii} = G;
    G = 0; 
end
%% Setting up small g and small b

smallG = YbusG;
for ii =1:3 
    smallG(ii,ii) = -smallG(ii,ii);
end
smallG = -smallG;

smallB = YbusB;
for ii =1:3 
    smallB(ii,ii) = -smallB(ii,ii);
end
smallB = -smallB;
%% Real Power Line Flow
Pmatrix = zeros(3,3);
for kk = 1:3
    for ll = 1:3
    P = V(kk)^2*(smallG(kk,ll)) - V(kk)*V(ll)*(smallG(kk,ll)*cos(D(ll)-D(kk)) - smallB(kk,ll)*sin(D(ll)-D(kk)));
    Pmat(kk,ll) = P; 
    end
end
%% Reactive Power Line Flow
Qmatrix = zeros(3,3);
bshunt = [0,.0775,.1275;.0775,0,.0775;.0775,.1275,0];
for kk = 1:3
    for ll = 1:3
        Q = -V(kk)^2*(smallB(kk,ll)+bshunt(kk,ll)/2) + V(kk)*V(ll)*(-smallG(kk,ll)*sin(D(kk) - D(ll)) + smallB(kk,ll)*cos(D(kk)-D(ll)));
        Qmat(kk,ll) = Q;
    end 
end
%% Calculating Sij for the Thermal Limits

S12 = sqrt(Pmat(1,2)^2+Qmat(1,2)^2);
S13 = sqrt(Pmat(1,3)^2+Qmat(1,3)^2);
S23 = sqrt(Pmat(2,3)^2+Qmat(2,3)^2);
S21 = sqrt(Pmat(2,1)^2+Qmat(2,1)^2);
S31 = sqrt(Pmat(3,1)^2+Qmat(3,1)^2);
S32 = sqrt(Pmat(3,2)^2+Qmat(3,2)^2);
%% Setting up output for fmincon


ceq = [Pvec{1}-x(1)/100;Pvec{2}-x(2)/100;Pvec{3}+3;Qvec{1} - x(3)/100; Qvec{2} - x(4)/100;Qvec{3}+1];
c=[S12*100-250;S13*100-180;S23*100-250;S21*100-250;S31*100-180;S32*100-250];


end