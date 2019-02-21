%% Problem 1
clear, clc
a = 5;
p = 10;
w = 3;
space = linspace(0,4*a,1000);
Mz1 = @(x) .75*p*x+.5*w*a*x;
Mz2 = @(x) -.25*p*x+.5*w*a*x+p*a;
Mz3 = @(x) .75*p*x+.5*a*w*x+p*a+.5*w*x^2+2*w*a^2-w*2*a*x;
Fy1 = @(x) .75*p+.5*w*a;
Fy2 = @(x) -.25*p+.5*w*a;
Fy3 = @(x) .75*p+.5*w*a-p-w*(x-2*a);
for i=1:length(space)/4
    M(i) = Mz1(space(i));
    Fy(i) = Fy1(space(i));
end
for i=length(space)/4+1:length(space)/2
    M(i) = Mz2(space(i));
    Fy(i) = Fy2(space(i));
end
for i = length(space)/2+1:length(space)
     M(i) = Mz3(space(i));
     Fy(i) = Fy3(space(i));
end
hold on
plot(space,M)
plot(space,Fy)

%% Problem 2
clear all, clc

clear, clc
a = 5;
p = 10;
w = 3;
space = linspace(0,3*a,999);
Mz1 = @(x) .5*w*x^2;
Mz2 = @(x) p*a*.5-.75*w*a^2;
Mz3 = @(x) 2.5*p*a-.75*w*a^2;
Fy1 = @(x) -w*x;
Fy2 = @(x) .5*p+.25*w*a;
Fy3 = @(x) -.5*p+.25*w*a;
for i=1:length(space)/3
    M(i) = Mz1(space(i));
    Fy(i) = Fy1(space(i));
end
for i=length(space)/3+1:2*length(space)/3
    M(i) = Mz2(space(i));
    Fy(i) = Fy2(space(i));
end
for i = 2*length(space)/3+1:length(space)
     M(i) = Mz3(space(i));
     Fy(i) = Fy3(space(i));
end
hold on
plot(space,M)
plot(space,Fy)

%% Problem 3
clear all, clc

clear, clc
fo = 5;
L = 10;
space = linspace(0,L,1000);
Mz1 = @(x) (-1/8)*fo*L^2+.5*fo*x^2;
Mz2 = @(x) 0;

Fy1 = @(x) -fo*L/2+fo*x;
Fy2 = @(x) 0;

for i=1:length(space)/2
    M(i) = Mz1(space(i));
    Fy(i) = Fy1(space(i));
end
for i=length(space)/2+1:length(space)
    M(i) = Mz2(space(i));
    Fy(i) = Fy2(space(i));
end

hold on
plot(space,M)
plot(space,Fy)

%% Problem 4
clear all, clc

clear, clc
lo = 5;
L = 10;
l = @(x) x*lo*(1-(x/L)^2);

Vy = @(x) -2*lo*L/3+lo*x-lo*x^3/(3*L^2);
Vyc = @(c) -2*lo*L/3+lo*c*L-lo*c^3*L/3;
Mz = @(x) lo*L^2/4-2*x*lo*L/3+x^2*lo/2-lo*x^4/(12*L^2);
Mzc = @(c) lo*L^2/4-2*c*lo*L^2/3+c^2*L^2*lo/2-lo*c^4*L^2/12;
hold on
yyaxis right
ylabel('Internal Shear Force [N]')
fplot(Vyc,[0,1],'linewidth',1.75)
yyaxis left
ylabel('Internal Moment [Nm]')
fplot(Mzc,[0,1],'linewidth',1.75)
title('Shear Force and Moment Inside the "Wing"')
legend('Shear Force [N]','Moment [Nm]')