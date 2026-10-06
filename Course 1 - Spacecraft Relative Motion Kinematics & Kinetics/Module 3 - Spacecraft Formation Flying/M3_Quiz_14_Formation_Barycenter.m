%% 
clear all; close all; clc;
import astro.*   % COE2RV / RV2COE from lib/+astro (run `startup` if not found)
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Quiz setup:
global mu
mu = 3.986e14; % mu Earth

% Chief orbit:
a     = 7000e3;           % [m]
e     = 0.01;
i     = 78*pi/180;        % [rad]
RAAN  = 120*pi/180;       % [rad]
omega = 33*pi/180;       % [rad]
f     = 45*pi/180;     % [rad] arbitrary true anomaly

%% Q1: Center-of-mass Cartesian center:
[rc_N, vc_N] = COE2RV(a,e,i,RAAN,omega,f,mu)

% Deputy orbit:
delta_f = 1*pi/180; % assume the Deputy is the leader.
[rd_N, vd_N] = COE2RV(a,e,i,RAAN,omega,f + delta_f,mu)

rCCM_N = (rc_N + rd_N) / 2 % No mass provided and assuming same mass between two satellites.
vCCM_N = (vc_N + vd_N) / 2


vpa(rCCM_N)'/1000

% ANS:
% # adjust the return matrix values as needed
% def result():
%     rCCM_N = [-1919.1629, 492.0015, 6661.9588]  # km
% 
%     return rCCM_N
% 
% result()

%% Q2: OE Barycenter:
% Since no mass is provided, the average f is simply in between the chief and deputy:
f_bc = (f + (f + delta_f))/2;
r_bc_N = COE2RV(a,e,i,RAAN,omega,f_bc,mu)

vpa(r_bc_N)'/1000

% ANS:
% # adjust the return matrix values as needed
% def result():
%     rOECM_N = [-1919.2372, 492.0234, 6662.2010] # km
% 
%     return rOECM_N
% 
% result()

%% Q3: Find the orbital period absolute difference between the Cart. COM & OE Barycenter:
period_c_d_bc = 2*pi*sqrt(a^3/mu)

[a_COM,~,~,~,~,~] = RV2COE(rCCM_N,vCCM_N,mu)
period_COM = 2*pi*sqrt(a_COM^3/mu)

period_c_d_bc - period_COM % Ans already in seconds.
