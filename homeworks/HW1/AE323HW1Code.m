%% Lots of variables
%% Problem 1a
clear
clc
N = @(b) (3-2*cot(b))/(3+2*tan(b));
%% Problem 1b
B = fsolve(N,.5)
%% Problem 1c
etamat = linspace(0,.25,100);
for i=1:length(etamat)
    N =@(b) N(b)-.25/100;
    betac(i) = fsolve(N,.01)
end
%% Problem 1cp2
plot(etamat,betac,'linewidth',1.5)
legend('Beta Values')
title('connorl2-Beta vs. Eta')
xlabel('Eta')
ylabel('Beta')
%% Problem 2c
clear
clc
syms b q O real
Tn = 2*b*q*sin(O)*cos(O)+q*sin(O)^2;
dTndO = diff(Tn,O)
dt2dO2 = diff(dTndO,O)
%% Problem 2e
clear
clc
B = @(O) -sin(O)*cos(O)/(cos(O)^2-sin(O)^2);
Bmat = linspace(0,6,1000);
for i=1:length(Bmat)
    B = @(O) B(O)-6/1000;
    O(i) = (360/(2*pi))*fsolve(B,1);
end
%%
plot(Bmat,O,'linewidth',1.5)
legend('ThetaStar Values')
title('connorl2-ThetaStar vs. Beta')
xlabel('Beta')
ylabel('ThetaStar')