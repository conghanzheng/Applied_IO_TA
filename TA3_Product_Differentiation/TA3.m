%% Exercise 7.2, (1)-(2)
% Price competition

% Preliminaries
clear all; clc;
syms t c0 c1 p0 p1 pi x; 
assume(t>0); 
assume(c0>=0); assume(c1>=0);
assume(p0>=0);assume(p1>=0);
assume(pi>=0);

% Solve profit maximization

% Indifferent consumer
xtilde = solve(-p0-t*x == -p1-t*(1-x),x);

% firm 0
f0 = (p0-c0)*xtilde;
grad0 = gradient(f0,p0);
hess0 = gradient(grad0,p0);
br0 = solve(grad0==0, p0, ReturnConditions=true); % best response 0

% firm 1
f1 = (p1-c1)*(1-xtilde);
grad1 = gradient(f1,p1);
hess1 = gradient(grad1,p1);
br1 = solve(grad1==0, p1, ReturnConditions=true); % best response 1

% Equilibrium prices, by solving simultaneous equations
p = solve(p0 == br0.p0, p1 == br1.p1, [p0,p1], ReturnConditions=true);

% Update xtilde
xtilde_new = solve(-p.p0-t*x == -p.p1-t*(1-x),x);

% Profits
pi0 = simplify((p.p0-c0)*xtilde_new);
pi1 = simplify((p.p1-c1)*(1-xtilde_new));

fprintf("\n ==== Exercise 7.2 (1) ==== \n")
fprintf("The second order cross partial derivative ")
fprintf('When t satisfies that ')
disp(subs(solve(hess0 < 0, t, ReturnConditions=true).conditions,"x","t"))
if solve(hess1 < 0, t, ReturnConditions=true).conditions ~= solve(hess0 < 0, t, ReturnConditions=true).conditions
    fprintf('and')
    disp(subs(solve(hess1 < 0, t, ReturnConditions=true).conditions,'x','t'))
end
fprintf('The Bertrand-Nash equilibrium: \n')
fprintf('p0* = ')
disp(p.p0)
fprintf('p1* = ')
disp(p.p1)
fprintf('pi0 = ')
disp(pi0)
fprintf('pi1 = ')
disp(pi1) 

% Second-order cross partial derivatives
dpi0_0 = diff(pi0,c0);
d2pi0_01 = diff(dpi0_0,c1);

dpi0_1 = diff(pi0,c1);
d2pi0_10 = diff(dpi0_1,c0);

dpi1_0 = diff(pi1,c0);
d2pi1_01 = diff(dpi1_0,c1);

dpi1_1 = diff(pi1,c1);
d2pi1_10 = diff(dpi1_1,c0);

fprintf("\n ==== Exercise 7.2 (2) ==== \n")
fprintf("The second order cross partial derivatives: \n ")
if d2pi0_01 == d2pi0_10
    fprintf("d^2(pi0)/[d(c0)d(c1)] = ")
    disp(d2pi0_01)
else
    fprintf("Symmetry of second derivatives is violated (pi0), check on that.")
end

if d2pi1_01 == d2pi1_10
    fprintf("d^2(pi1)/[d(c0)d(c1)] = ")
    disp(d2pi1_01)
else
    fprintf("Symmetry of second derivatives is violated (pi1), check on that.")
end

%% Exercise 7.3

syms p_i t x p n c f
assume(p>=c)
assume(in(t,'real') & t>0)
assume(in(f,'real') & f>0)
assume(in(n,'real') & n>0)

fprintf("\n ==== Exercise 7.3 ==== \n")

% Indifferent consumer
xtilde = solve(p_i + t*(x^2) == p + t*((1/n-x)^2), x);
xtilde = simplify(xtilde);
fprintf("Indifferent consumer xtilde = ")
disp(xtilde)

syms xtilde(y,z)
xtilde(y,z) = solve(y + t*(x^2) == z + t*((1/n-x)^2), x);

% Demand 
syms D(y,z)
D(y,z) = 2*xtilde(y,z); 
fprintf("Demand D = ")
disp(D(p_i,p))

% Profit
pi = simplify((p_i - c)*D(p_i,p) -  f);
grad = diff(pi, p_i);
hess = diff(pi, p_i, p_i);
br = solve(grad==0, p_i, ReturnConditions=true);

p_symeq = solve(br.p_i == p, p, ReturnConditions=true);
fprintf("Symmetric equilibrium price p_i = p = ")
disp(p_symeq.p)

fprintf("The profic margin p-c =")
disp(simplify(p_symeq.p-c))

pi_symeq = simplify((p_symeq.p - c)*D(p_symeq.p,p_symeq.p) - f);
fprintf("Profit pi = ")
disp(pi_symeq)
n_free = solve(pi_symeq == 0, n, ReturnConditions=true).n; 
fprintf("Market size under free entry n^c = ")
disp(n_free)

% Social planner problem
syms d g(x)
g(x) = 1/(1/(2*n)-0); % uniform distribution 1/(b-a)

% E(d^2)
E = int(x^2*g(x), x, 0, 1/(2*n));
fprintf("E(d^2) = ")
disp(E)

% The social planner wants to minimise the unproductive expenses (fixed and transport costs)
cost = n*f + t*E;
grad = diff(cost, n);
hess = diff(cost, n, n);
n_plan = solve(grad==0, n, ReturnConditions=true).n;

fprintf("The optimal amount of firms n^* = ")
disp(n_plan)

if isAlways(n_free > n_plan)
    fprintf("and n^* < n^c \n")
else 
    fprintf("Error, please check.")
end