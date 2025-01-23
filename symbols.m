
clear all; clc;
syms theta real
syms omega real
syms u real
syms b real
syms k real
syms J
x = [theta; omega];
    
f1 = omega;
f2 = (u - b*omega - k*theta)/J;
f = [f1; f2];
h = theta;
E = solve([theta == pi * 17/36, f1 == 0, f2 == 0], 'Real', true);

A = jacobian(f,x);
B = jacobian(f,u);
C = jacobian(h,x);
D = jacobian(h,u);

psi = [-1, 2, 1.5, -2];
L0 = 2;
b = .2;
k = 20;
u = double(subs(E.u));
theta = double(E.theta);
omega = double(E.omega);
L = [0.5, 0.01, 0.01, 0.08];
J = L0 + sum(L .* cos((1:4) * theta + psi));

A = double(subs(A));
B = double(subs(B));
C = double(subs(C));
D = double(subs(D));

s = tf('s');
GG=C*(inv(s*eye(2)-A))*B+D

syms xi_star real
syms S real
xi = solve(S == exp(-pi*xi_star/sqrt(1-xi_star^2)), xi_star, 'Real', true);
S = 0.08;
xi = double(subs(xi))

%% regolatore dinamico
syms G_e
syms M_star real
M = solve(abs(G_e) + 20*log10(M_star)==0,M_star,'Real',true)
