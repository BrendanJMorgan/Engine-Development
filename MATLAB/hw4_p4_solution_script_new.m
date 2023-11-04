% MAE 150/250P HW 4 P4
% Fall 2023
% Units in: SI
clear all;
close all;
clc;
% constants
R = 287;        % J/kg-K, gas constant for Air
g = 9.81;       % m/s^2
gamma = 1.4;    

M0 = 0.65;
T0 = 230;
p0 = 35e3;
v0 = M0*sqrt(gamma*R*T0);
cp = 1004;
alpha_opt = .91;
eta_prop = 0.78;
eta_gb = 0.99;

pi_d = 0.94;
pi_c = 30;
e_c = 0.95;
mdot_o = 33;

Tt4 = 1500;
pi_b = 0.98;
eta_b = 0.98;
Qr = 43800e3;

e_hpt = 0.82;
eta_mhpt = 0.97;
eta_mlpt = 0.97;
eta_lpt = .92;
eta_t = 0.9;

eta_n = 0.99;
p9 = p0;
%% find
Tt0 = T0*(1+(gamma-1)/2*M0^2);
pt0 = p0*(1+(gamma-1)/2*M0^2)^(gamma/(gamma-1));
Tt2 = Tt0;
pt2 = pt0*pi_d;
pt3 = pt2*pi_c;
tc = pi_c^((gamma-1)/(gamma*e_c));
Tt3 = Tt2*tc;
pt4 = pt3*pi_b;
f = (Tt4-Tt3)/(Qr*eta_b/cp-Tt4);
Tt45 = Tt4-(Tt3-Tt2)/(eta_mhpt*(1+f));
t_hpt = Tt45/Tt4;
pi_hpt = t_hpt^(gamma/(gamma-1)/e_hpt);
pt45 = pt4*pi_hpt;
Power_lpt = mdot_o*(1+f)*eta_lpt*alpha_opt*cp*Tt45*(1-(p9/pt45)^((gamma-1)/gamma));
Fprop = eta_prop*eta_gb*eta_mlpt*Power_lpt/v0;

% Tt5 = Tt45-Power_lpt/mdot_o/cp;
% t_lpt = Tt5/Tt45;
% 
% syms EE
% e_lpt = double(vpasolve(eta_lpt == (1-t_lpt)/(1-t_lpt^(1/EE)),EE));
% 
% pi_lpt = t_lpt^(gamma/(gamma-1)/e_lpt);
% pt5 = pt45*pi_lpt;


V9 = sqrt(2*(1-alpha_opt)*eta_n*cp*Tt45*(1-(p9/pt45)^((gamma-1)/gamma)));
Fcore = mdot_o*((1+f)*V9-v0);

PSFC = mdot_o*f/(eta_gb*eta_mlpt*Power_lpt+(mdot_o*.5*((1+f)*V9^2-v0^2)));
eta_th = 1/(f*mdot_o*Qr)*(eta_gb*eta_mlpt*Power_lpt+(mdot_o*.5*((1+f)*V9^2-v0^2)));
eta_p = (Fcore+Fprop)*v0/(eta_gb*eta_mlpt*Power_lpt+(mdot_o*.5*((1+f)*V9^2-v0^2)));
eta_tot = eta_th*eta_p;
KE_core = (mdot_o*.5*((1+f)*V9^2-v0^2));
P_prop = eta_gb*eta_mlpt*Power_lpt;
%PSFC = mdot_o*f/(P_prop+KE_core);