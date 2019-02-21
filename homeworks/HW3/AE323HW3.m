%% AE323 HW3 P1
clear
clc

t = .005;
g = 10;
p = 7800;
b = 12*t;
a = 6*t;
x = 0;
z = .5*b+t;
sigmaxx = @(L) z*p*g*(t*(2*a+b))*(.5*L^2-L*x+.5*x^2)/(b^3*t/12+b^2*t*a/2+t^2*b*a)
fsolve(@(L) sigmaxx(L)-2*10^8,10)
fplot(sigmaxx,[0,10])
%% 
clear
clc
syms p g b t a x L
z = -(t+.5*b)
sigmaxxdes = simplify(12*(t+.5*b)*p*g*(2*a+b)*(-.5*x^2+x*L-L^2/2)/(b^3+6*a*b^2))
sigmaxxcon = simplify(z*p*g*(t*(2*a+b))*(.5*L^2-L*x+.5*x^2)/(b^3*t/12+b^2*t*a/2+t^2*b*a))
simplify(sigmaxxcon/sigmaxxdes)

%% Problem 2
clear
clc

syms t
a = 12*t;
b = 10*t;
At = t*a;
Ab = b*t;
yt = 0;
yb = -t/2-b/2;
Astar = At+Ab;

neut = (At*yt+Ab*yb)/Astar

%% 
clear 
clc

syms a b t p L x
At  =t*a;
Ab = b*t;
yt = 0;
yb = -t/2-b/2;
Astar = At+Ab;
neutral = (At*yt+Ab*yb)/Astar
A1cold = 0;
A2cold = -t/2-b/2;
A1cnew = 0-neutral;
A2cnew = A2cold-neutral;
Iyy = simplify((A1cnew)^2*At+a*t^3/12+(A2cnew)^2*Ab+b^3*t/12);
A1cnew = subs(A1cnew,{a,b},{12*t,10*t});
A2cnew = subs(A2cnew,{a,b},{12*t,10*t});
Iyy = simplify((A1cnew)^2*At+a*t^3/12+(A2cnew)^2*Ab+b^3*t/12);
Iyy = subs(Iyy,{a,b},{12*t,10*t});
oldtopz = t/2;
oldbotz = -.5*t-b;
newtopz = oldtopz-neutral
newtopz = subs(newtopz,{a,b},{12*t,10*t})
newbotz = oldbotz-neutral
newbotz = subs(newbotz,{a,b},{12*t,10*t})
% 
% My = p*a*(x*L-.5*L^2-.5*x^2);
% 
% sigxxtop = simplify(-My*newtopz/Iyy)
% sigxxbot = simplify(-My*newbotz/Iyy)
% 
% sigxxtop = subs(sigxxtop,{t,a,b,p,L},{.005,12*.005,10*.005,500000,10})
% sigxxtop
% fplot(@(x) sigxxtop(x))