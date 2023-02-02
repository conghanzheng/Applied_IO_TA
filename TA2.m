% Preliminaries
clear all; clc;
syms q1 q2 c c1; 
assume(c1>0); assume(c1<2*c);
assume(q1>=0); assume(q2>=0);
assume(q1<1); assume(q2<1); % because the market demand is q = 1-p

% Solve profit maximization

% firm 1
f1 = (1-q1-q2-c1)*q1;
grad1 = gradient(f1,q1);
hess1 = gradient(grad1,q1);
br1 = solve(grad1==0,q1, ReturnConditions=true); % best response: q1 as a function of q2 and c1

% firm 2
f2 = (1-q1-q2-(2*c-c1))*q2;
grad2 = gradient(f2,q2);
hess2 = gradient(grad2,q2);
br2 = solve(grad2==0,q2, ReturnConditions=true); % best response: q2 as a function of q1 and c2

% Solve for quantities, simultaneous equations
q = solve(q1 == br1.q1, q2 == br2.q2, ReturnConditions=true);

% Total profit
pi1 = simplify((1-q.q1-q.q2-c1)*q.q1);
pi2 = simplify((1-q.q1-q.q2-(2*c-c1))*q.q2);
pi = simplify(pi1+pi2);

% Total profit on c1
gradpi = gradient(pi,c1);
hesspi = gradient(gradpi,c1);
c1_solvepi = solve(gradpi==0,c1, ReturnConditions=true);

% Print results
if hess1 < 0 && hess2 <0 
    fprintf('q1* = ')
    disp(q.q1)
    fprintf('q2* = ')
    disp(q.q2)
    fprintf('Total profit = ')
    disp(pi)
else
    disp("No solution found for the two profix maximization problems.")
end

if hesspi > 0 
    fprintf("The total profit is convex in c1 and reaches its local minimum at c1 = ")
    disp(c1_solvepi.c1)
elseif hesspi < 0
    fprintf("The total profit is concave in c1 and reaches its local maximum at c1 = ")
    disp(c1_solvepi.c1)
else
    fprintf("The profit function doesn't necessarily have an extreme at c1, it's probably an inflection point, check on that.")
end