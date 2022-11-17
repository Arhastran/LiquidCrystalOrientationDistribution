%Simulations for LQ RI distribution using finite element method 
%Simple Laplacian method with angles, fields and elastic constants approimation



syms x y 
syms phi(x,y)
syms phi1(x)
syms phi2(y)
lapphi = laplacian(phi, [x y]) == 0
%lapphi = diff(phi1, x, 2) + diff(phi2,y,2) == 0;
%lapphi = diff(phi1,x,2)



cond1 = phi(0,y) == pi/2;
cond2 = phi(1,y) == pi/2;
cond3 = phi(x,0) == 0;
cond4 = phi(x,1) == 0;


conds = [cond1 cond2 cond3 cond4];

x1 = 0:0.1:1;
y1 = 0:0.1:1;
initial = ones(length(x1))*pi/2;



solphi = pdepe(2,lapphi,initial,@bcfun,x1,y1)

%V = odeToVectorField(lapphi)

solphi = dsolve(lapphi,conds,'IgnoreAnalyticConstraints',true);







