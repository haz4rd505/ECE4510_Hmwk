% function for injected Power for 3 bus system of economic dispatch code

function [c,ceq] = injectedpower(x)

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
    Pvec{ii} = F
    F=0;
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
    G=0;
end
%% Setting up output for fmincon

ceq = [Pvec{1}-x(1)/100;Pvec{2}-x(2)/100;Pvec{3}+3;Qvec{1} - x(3)/100; Qvec{2} - x(4)/100;Qvec{3}+1];
c=[];


end