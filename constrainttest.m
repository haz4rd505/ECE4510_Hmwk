

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

% x = case3
% x = case4
x = case5



V = [x(5), x(6),x(7)];
D = [0,x(8),x(9)];

Pvec={};
F=0;
for ii = 1:3
    for jj = 1:3
        P = V(ii)*V(jj)*(YbusG(ii,jj)*cos(D(ii) - D(jj)) + YbusB(ii,jj)*sin(D(ii)-D(jj)));
        F = F + P ;
    end
    Pvec{ii} = F;
    F = 0;
end

Qvec = {};
G=0;
for ii = 1:3
    for jj = 1:3
        Q = V(ii)*V(jj)*(YbusG(ii,jj)*sin(D(ii) - D(jj)) - YbusB(ii,jj)*cos(D(ii)-D(jj)));
        G = G + Q ;
        Q;
    end
    Qvec{ii} = G;
    G=0;
end


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


Pmatrix = zeros(3,3);
for kk = 1:3
    
    for ll = 1:3
    P = V(kk)^2*(smallG(kk,ll)) - V(kk)*V(ll)*(smallG(kk,ll)*cos(D(ll)-D(kk)) - smallB(kk,ll)*sin(D(ll)-D(kk)));
    Pmat(kk,ll) = P; 
    end
end
 
Qmatrix = zeros(3,3);
for kk = 1:3
    for ll = 1:3
        Q = -V(kk)^2*(smallB(kk,ll)) + V(kk)*V(ll)*(-smallG(kk,ll)*sin(D(kk) - D(ll)) + smallB(kk,ll)*cos(D(kk)-D(ll)));
        Qmat(kk,ll) = Q;
    end 
end

