figure(3)
fplot(@(x) (2*sqrt(3) + 2*sqrt(3)*cos(240*pi*x + 2*pi/3) + 2*sin(240*pi*x+2*pi/3)), [0 .05],'r')
hold on 
fplot(@(x) (2*cos(120*pi*x + pi/6)), [0 .05], 'k')
fplot(@(x) (4*cos(120*pi*x + pi/3)), [0 .05], 'b')
legend('p(t)','i(t)','v(t)')
title('Problem 4 Waveplot')
xlabel('Time (t)')
ylabel('Units in W/A/V')
%fplot(@(x) (2*sqrt(3)))

%%

fplot(@(x) (5*.75)^x/factorial(x)*exp(-(5*.75)))



%%