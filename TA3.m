% Preliminaries
clear all; clc;
syms t c0 c1 p0 p1 pi; 
assume(t>0); 
assume(c0>=0); assume(c1>=0);
assume(p0>=0);assume(p1>=0);
assume(pi>=0);

% Solve profit maximization

% firm 0
f0 = (p0-c0)*(p1-p0+t)/(2*t);
grad0 = gradient(f0,p0);
hess0 = gradient(grad0,p0);
br0 = solve(grad0==0, p0, ReturnConditions=true); % best response 0

% firm 1
f1 = (p1-c1)*(p0-p1+t)/(2*t);
grad1 = gradient(f1,p1);
hess1 = gradient(grad1,p1);
br1 = solve(grad1==0, p1, ReturnConditions=true); % best response 1

% Equilibrium prices, by solving simultaneous equations
p = solve(p0 == br0.p0, p1 == br1.p1, [p0,p1], ReturnConditions=true);

% Profits
pi0 = simplify((p.p0-c0)*(p.p1-p.p0+t)/(2*t));
pi1 = simplify((p.p1-c1)*(p.p0-p.p1+t)/(2*t));

% Second-order cross partial derivatives
dpi0_0 = diff(pi0,c0);
d2pi0_01 = diff(dpi0_0,c1);

dpi0_1 = diff(pi0,c1);
d2pi0_10 = diff(dpi0_1,c0);

dpi1_0 = diff(pi1,c0);
d2pi1_01 = diff(dpi1_0,c1);

dpi1_1 = diff(pi1,c1);
d2pi1_10 = diff(dpi1_1,c0);

% Print results
fprintf('When t satisfies that ')
disp(solve(hess0 < 0, t, ReturnConditions=true).conditions)
if solve(hess1 < 0, t, ReturnConditions=true).conditions ~= solve(hess0 < 0, t, ReturnConditions=true).conditions
    fprintf('and')
    disp(subs(solve(hess1 < 0, t, ReturnConditions=true).conditions,'x','t'))
end
fprintf('The Bertrand-Nash equilicrium prices: \n')
fprintf('Equilibrium price p0* = ')
disp(p.p0)
fprintf('Equilibrium price p1* = ')
disp(p.p1)
fprintf('Equilibrium profit pi0 = ')
disp(pi0)
fprintf('Equilibrium profit pi1 = ')
disp(pi1)  

if d2pi0_01 == d2pi0_10
    fprintf("The second order cross partial derivative ")
    fprintf("d^2(pi0)/[d(c0)d(c1)] = ")
    disp(d2pi0_01)
else
    fprintf("Symmetry of second derivatives is violated (pi0), check on that.")
end

if d2pi1_01 == d2pi1_10
    fprintf("The second order cross partial derivative ")
    fprintf("d^2(pi1)/[d(c0)d(c1)] = ")
    disp(d2pi1_01)
else
    fprintf("Symmetry of second derivatives is violated (pi0), check on that.")
end