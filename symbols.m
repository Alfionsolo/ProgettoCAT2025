clear; clc;
syms theta real
syms omega real
syms u real
syms b real
syms k real
syms j real
syms T(theta) real
%syms alpha real
x = [theta; omega];

f1 = omega;
f2 = (u*T - b*omega - k*theta)/j;
f = [f1; f2];
h = theta;
E = solve([theta == pi * 5/9,f1 == 0,f2 == 0], 'Real' , true);

%%
A = jacobian(f,x);
B = jacobian(f,u);
C = jacobian(h,x);
D = jacobian(h,u);

j=400;
b = 0.5;
k = 500;
alpha=pi*5/18;
theta = double(E.theta);
omega = double(E.omega);

T = cos(alpha)/(1-(sin(alpha)*cos(theta))^2);

dT=-(2*cos(alpha)*sin(alpha)^2*cos(theta)*sin(theta))/(1-sin(alpha)^2*cos(theta)^2)^2;

u = double(subs(E.u)); %In B, u deve essere un simbolo
df2x1Num=-(k-u* dT )/j;
%% Linearizzazione
A = double(subs(A));
B = double(subs(B));
C = double(subs(C));
D = double(subs(D));

s = tf('s');
GG=C*(inv(s*eye(2)-A))*B+D;

P=[1,0.00125,1.25];
y_val=roots(P);
%-0.000625000000000 + 1.118033814057071i
%%


%%PROVE E CONTRO PROVE
% %E = solve([theta == pi * 5/9,f1 == 0,f2 == 0], 'Real' , true);
% %T = cos(alpha)/(1-(sin(alpha)*cos(theta))^2);
% 
% %
% [df1x1]=gradient(f1,theta); 
% [df1x2]=gradient(f1,omega); 
% [df2x1]=gradient(f2,theta); 
% [df2x2]=gradient(f2,omega); 
% 
% E = solve([theta == pi*5/9,f1 == 0,f2 == 0],'Real',true);
% 
% %T = cos(alpha)/(1-(sin(alpha)*cos(theta))^2);
% d=diff(T,theta);
% %dT=-(2*cos(alpha)*sin(alpha)^2*cos(theta)*sin(theta))/(sin(alpha)^2*cos(theta)^2 - 1)^2;
% 
% j=400;
% b = 0.5;
% k = 500;
% alpha=pi*5/18;
% theta = double(E.theta);
% omega = double(E.omega);
% theta=pi * (5 / 9); %1OO
% %u= pi* 424.5; %kx1e/tao
% 
% T = cos(alpha)/(1-(sin(alpha)*cos(theta))^2);
% dT=-(2*cos(alpha)*sin(alpha)^2*cos(theta)*sin(theta))/(sin(alpha)^2*cos(theta)^2 - 1)^2;%La puttana
% 
% u = double(subs(E.u)); %In B, u deve essere un simbolo
% 
% df2x1Num=-(k-u* ddd )/j;
% A = [df1x1,df1x2; df2x1Num, df2x2];
