%ECE 4510 Economic Dispatch Project 

%These are our constraints for case 1
Vmin = .94;
Vmax = 1.06;
Pmin1=0;
Pmin2=0;
Pmax1 = 200;
Pmax2 = 270;
Qmin1=-100;
Qmin2=-100;
Qmax1=100;
Qmax2=100;
%% Case 1 Minimizing Production Cost with a lower bound on the Real Power
% Generation
Aeq = [1,1]; beq = 300;
lb = [0 0]
 x0 = [0 0]
fmincon(@(x) .1225*x(1)^2 + x(1) + .81*x(2)^2 + 5*x(2) + 335 + 50,x0,[],[],Aeq,beq,lb)



%% Case 2 Minimizing Production Cost with an upper and lower bound on the
% Real Power Generation
Aeq = [1,1]; beq = 300;
lb = [0 0 ]
ub = [200 270];
 x0 = [0 0];
fmincon(@(x).1225*x(1)^2 + x(1) + .81*x(2)^2 + 5*x(2) + 335 + 50,x0,[],[],Aeq,beq,lb,ub)   

%% Case 3 Minimizing Production Costs while using the Power flow equations
% as constraints and lower and upper bound for real power generation.
% Ybus matrix includes the shunts

x0 = [ 0 0 0 0 1 1 1 0 0];
lb = [0 0 -100 -100 0 0 0 -Inf -Inf]
ub = [200 270 100 100 Inf Inf Inf Inf Inf]
nonlcon = @injectedpower
case3 = fmincon(@(x).1225*x(1)^2 + x(1) + .81*x(2)^2 + 5*x(2) + 335 + 50,x0,[],[],[],[],lb,ub,nonlcon) 
case3out = case3;
case3out(8) = case3(8)*360/(2*pi)
case3out(9) = case3(9)*360/(2*pi)
%% Case 4 Minimizing Production with Power flow constraints, Thermit Limit Constraints, and voltage constraints.


Aeq = []; beq = []
x0 = [ 0 0 0 0 1 1 1 0 0];
lb = [0 0 -100 -100 .94 .94 .94 -Inf -Inf]
ub = [200 270 100 100 1.06 1.06 1.06 Inf Inf]
nonlcon = @allpowerconstraints
case4 = fmincon(@(x).1225*x(1)^2 + x(1) + .81*x(2)^2 + 5*x(2) + 335 + 50,x0,[],[],[],[],lb,ub,nonlcon)
case4out = case4;
case4out(8) = case4(8)*360/(2*pi)
case4out(9) = case4(9)*360/(2*pi)

%% Case 5, same as Case 4 with an 8% increase in load
% Load increase is included in pconsextraload


Aeq = []; beq = []
x0 = [ 0 0 0 0 1 1 1 0 0];
lb = [0 0 -100 -100 .94 .94 .94 -Inf -Inf]
ub = [200 270 100 100 1.06 1.06 1.06 Inf Inf]
nonlcon = @pconsextraload
case5 = fmincon(@(x).1225*x(1)^2 + x(1) + .81*x(2)^2 + 5*x(2) + 335 + 50,x0,[],[],[],[],lb,ub,nonlcon)
case5out = case5;
case5out(8) = case5(8)*360/(2*pi)
case5out(9) = case5(9)*360/(2*pi)














