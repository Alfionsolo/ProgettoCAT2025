
%% Progetto Controlli Automatici T B3
%GRUPPO 15- Gaia Margherita, Taruffi Alice,  Finocchiaro Alfio.
clear all; clc;

%% Parametri dati
th_e = pi * (5 / 9); %1OO
beta = 0.5;
alpha= pi* (5/18);
x_e = [th_e; 0];
u_e = pi* 424.5; %kx1e/tao
k = 500;
J = 400;
T = @(theta) cos(alpha)/(1-(sin(alpha)*cos(theta))^2);

%% Punto 1 - Linearizzazione
dT= @(theta)-(2*cos(alpha)*sin(alpha)^2*cos(theta)*sin(theta))/(1-sin(alpha)^2*cos(theta)^2 )^2;
%OPPURE:
%dT= @(th) -(cos(alpha)*sin(alpha)^2*sin(2*th)) / (sin(alpha)^2 *cos(th)^2-1)^2;

%Valutata:
%dT= -(2*cos(alpha)*sin(alpha)^2*cos(th_e)*sin(th_e))/(1-sin(alpha)^2*cos(th_e)^2 )^2;

%COMMENTARE PER VALUTARE dT all'equilibrio
df2dx1 = (u_e*dT(x_e(1))/J)-(k/J);
%SCOMMENTARE PER VALUTARE dT in th_e:
%df2dx1 = (u_e*dT/J)-(k/J);

%che sarebbe:
%df2x1=-(k-u* dtao )/j
%df2dx1= - k/400 -(u*cos(alpha)*sin(alpha)^2*cos(theta)*sin(theta))/(200*(sin(alpha)^2*cos(theta)^2 - 1)^2);  

df2dx2 = - beta/J;
df2du=T(x_e(1))/J;

%linearizzazione nell'intorno di equilibrio
A = [0, 1; df2dx1, df2dx2];
B = [0; df2du];
C = [1 0];
D = 0;

%% Punto2 - Funzione di trasferimento
s=tf('s');
G=C*(inv(s*eye(2)-A))*B+D;

%% Patch zone proibite
patched(2, G);

%% Regolatore statico
% mu_s = R(0) = L(0)/G(0)
% Per la specifica 5, L(0) > 10 ^ 50/20
% (conversione da dB a valore assoluto)
G0 = abs(evalfr(G,0));
L0 = 10^(50/20);

% Sistema esteso
% mu_s > L0/G0 ~ 6300;
mu_s = 7000;
G_e = mu_s * G;

%% Regolatore dinamico
% scenario B => rete anticipatrice
omega_c_star = 330;
M_f_star = 75;

% la funzione bode dà modulo e argomento di G_e(j omega_c*)
[mag_G_star, arg_omega_G_star, omega_c_star] = bode(G_e, omega_c_star);

% M* e phi*
M_star = 1/mag_G_star;
phi_star = deg2rad(M_f_star - 180 - arg_omega_G_star);
if M_star <= 1
    M_star %#ok<NOPTS>
    return
end
if phi_star < 0 || phi_star > pi/2
    phi_star %#ok<NOPTS>
    return
end
% formule di inversione
alpha_tau = (cos(phi_star)-1/M_star)/(omega_c_star*sin(phi_star));
tau = (M_star - cos(phi_star))/(omega_c_star*sin(phi_star));
if cos(phi_star) < 1/M_star
    disp('alpha negativo')
    return
end
% fine rete anticipatrice
R_a = (1+tau*s)/(1+alpha_tau*s);

% Per rispettare la specifica sul disturbo di misura inseriamo un polo alle
% alte frequenze
R_hf = 1 / (1 + s/0.8e4);

% Regolatore dinamico completo
R_d = R_a * R_hf;
L = G_e * R_d;

%% Diagrammi di Bode di Ge con specifiche
patched(3, G_e);
patched(4, L);

%% Punto 4 - Test sistema linearizzato
% 1 - Risposta al gradino
F = L / (1+L);
W = 3;
T_sim = 2;
[y_step,t_step] = step(W*F, T_sim);
W_max = W * (1+0.08);
T_star = 0.1;
figure(6)
hold on
plot(t_step,y_step,'b');
% vincolo sovraelongazione
patch([0,T_sim,T_sim,0],[W_max,W_max,W+1,W+1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
% vincolo tempo assestamento
patch([T_star,T_sim,T_sim,T_star],[W*(1-0.01),W*(1-0.01),0,0],'g','FaceAlpha',0.3,'EdgeAlpha',0.5);
patch([T_star,T_sim,T_sim,T_star],[W*(1+0.01),W*(1+0.01),W+1,W+1],'g','FaceAlpha',0.3,'EdgeAlpha',0.5);
legend(["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"]);
ylim([W-1,W+1]);
hold off
% 2 - Risposta a d(t)
S = 1 / (1+L);
d = @(t) 3*sum(sin((1e-2)*(1:3).*t))/10^(50/20);
tt = (0:1e-2:1e3);
dd = arrayfun(d,tt);
y_d = lsim(S,dd,tt);
figure(8)
hold on
plot(tt,dd,'m')
plot(tt,y_d,'b')
legend('dd','y_d')
hold off
% 3 - Risposta a n(t)
n = @(t) 3*sum(sin((1e+5)*(1:3).*t))/10^(60/20);
tt = (0:1e-5:1);
nn = arrayfun(n,tt);
y_n = lsim(-F,nn,tt);
figure(9)
hold on
plot(tt,nn,'m')
plot(tt,y_n,'b')
legend('nn','y_n')
xlim([0,1e-3])
hold off

y_d_test = lsim(S,dd,tt);
y_n_test = lsim(-F,nn,tt);
y_w_test = 3*step(F,tt);
y = y_w_test + y_d_test + y_n_test;
figure(10)
hold on
rif = 3 * heaviside(tt);
plot(tt,rif,'b');
plot(tt,y,'r');
legend('w','y');
xlim([0,0.3])
ylim([W-1,W+1])
hold off
%% Punto 3 - Mappatura specifiche
function patched(n, tf)

% 1) Errore a regime |e∞| ≤ e⋆ = 0.05 in risposta a un gradino w(t) = 3(t) e d(t) = 3(t)
e_star = 0.05;
W = 3;
d = 3;
mu_s = (W + d) / e_star;


% 2) Margine di fase >= 45°
Mf_spec = 45;

% 3) Il sistema può accettare una sovraelongazione percentuale al massimo dell’8% : S% ≤ 8%
S_max = 0.08;
logsq = (log(S_max))^2;
xi = sqrt(logsq/(pi^2+logsq));
Mf_spec = max(Mf_spec, xi*100);

% 4) Il tempo di assestamento alla ε% = 1% deve essere inferiore al valore fissato: Ta,ε = 0.1s
T_star = 0.1;
epsilon = 0.01;
omega_c = 460.5 / (T_star * Mf_spec);
omega_Ta_MAX=omega_c;

% 5) Il disturbo sull' uscita d(t), con una banda limitata nel range di pulsazioni [0, 0.1], deve essere abbattuto di almeno 50 dB.
A_d = 50;
omega_d_min = 1e-4;
omega_d_MAX = 0.1;

% 6) Il rumore di misura n(t), con una banda limitata nel range di pulsazioni [1.5·10^4 , 10^6], deve essere abbattuto di almeno 60 dB.
A_n = 60;
omega_n_min = 1.5e4;
omega_n_MAX = 1e6;

% Dati plot
omega_plot_min = 1e-4;
omega_plot_max = 1e6;
xs = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
ys = [A_d; A_d; -150; -150];

%Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot
omega_Ta_min = 1e-4; % lower bound per il plot

Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -300; -300];

% Specifiche su n (massima pulsazione di taglio)
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];


figure(n);
hold on
patch(xs, ys,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
%Legenda colori
Legend_mag = ["A_d", "A_n", "\omega_{c,min}"];
legend(Legend_mag);

margin(tf, {omega_plot_min, omega_plot_max});
grid on; zoom on;
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

%Legenda colori
Legend_arg = ["Ge(j\omega)"; "M_f"];
legend(Legend_arg);
hold off
end



