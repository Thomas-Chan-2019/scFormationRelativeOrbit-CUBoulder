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
a     = 10000e3;           % [m]
e     = 0.2;
i     = 37*pi/180;        % [rad]
RAAN  = 40*pi/180;       % [rad]
omega = 65*pi/180;       % [rad]
f0    = 10*pi/180;     % [rad] arbitrary true anomaly
f1    = f0;
f2    = f0 + 60*pi/180;

theta_0 = omega + f0;
theta_2 = omega + f2;

del_a = 0;
del_e = 0.0001;
del_i = 0.001; % [rad]
del_RAAN = 0.001; % [rad]
del_omega = 0.001; % [rad]
del_M0 = -0.001; % [rad]

%% Quiz 16 - Q1: OE Difference Solution -> emall e:
eta = sqrt(1 - e^2); 

% Use Mean Anomaly Drift Evaluation (M3_Quiz_14, only in notes!):
M0 = E2M(f2E(f0,e),e); M1 = M0;
M2 = E2M(f2E(f2,e),e);
del_M2 = del_M0 - 3/2 * (M2 - M0) * del_a/a;

% r = a*(1 - e*cos(f2));

del_x = sqrt((e*del_M2/eta)^2 + (del_e + del_a/a)^2);
del_y = sqrt(4*del_e^2 + e^2*(del_M2/eta - del_omega - cos(i)*del_RAAN)^2);
del_z = sqrt(del_i^2 + sin(i)^2 * del_RAAN^2);

f_x = atan2(e*del_M2, -eta*(del_e + del_a/a));
f_y = atan2(e*(del_M2/eta - del_omega - cos(i)*del_RAAN), -2*del_e);
f_z = atan2((cos(omega)*del_i + sin(omega)*sin(i)*del_RAAN),(sin(omega)*del_i - cos(omega)*sin(i)*del_RAAN));
theta_z = atan2(del_i, -sin(i)*del_RAAN);

% rho_t_H:
x_1 = del_a + a*del_x*cos(f1 - f_x);
y_1 = a*(del_M2/eta + del_omega + cos(i)*del_RAAN) - a*del_y*sin(f1 - f_y) - a*e*sin(2*f1)*del_e/2;
z_1 = a*del_z*cos(theta_0 - theta_z) - a*e*del_z*sin(2*f1 - f_z)/2 - a*e*(sin(omega)*del_i - cos(omega)*sin(i)*del_RAAN)/2;

rho_1_H = vpa([x_1, y_1, z_1])/1000

% rho_t_H:
x_2 = del_a + a*del_x*cos(f2 - f_x);
y_2 = a*(del_M2/eta + del_omega + cos(i)*del_RAAN) - a*del_y*sin(f2 - f_y) - a*e*sin(2*f2)*del_e/2;
z_2 = a*del_z*cos(theta_2 - theta_z) - a*e*del_z*sin(2*f2 - f_z)/2 - a*e*(sin(omega)*del_i - cos(omega)*sin(i)*del_RAAN)/2;

rho_2_H = vpa([x_2, y_2, z_2])/1000

%% Solution on top is WRONG becuase e is much larger than assumed, pursuing another solution:
r1 = a*eta^2/(1 + e*cos(f1));
r2 = a*eta^2/(1 + e*cos(f2));

del_M1 = del_M0 - 3/2* (M1 - M0) * del_a/a

x1 = r1/a*del_a + a*e*sin(f1)/eta * del_M1 - a*cos(f1)*del_e;
y1 = r1/eta^3 * (1 + e*cos(f1))^2 * del_M1 + r1*del_omega + r1*sin(f1)/eta^2 * (2 + e*cos(f1))*del_e + r1*cos(i)*del_RAAN;
z1 = r1*(sin(theta_0)*del_i - cos(theta_0)*sin(i)*del_RAAN);

rho_1H_new = vpa([x1,y1,z1])/1000

x2 = r2/a*del_a + a*e*sin(f2)/eta * del_M2 - a*cos(f2)*del_e;
y2 = r2/eta^3 * (1 + e*cos(f2))^2 * del_M2 + r2*del_omega + r2*sin(f2)/eta^2 * (2 + e*cos(f2))*del_e + r2*cos(i)*del_RAAN;
z2 = r2*(sin(theta_2)*del_i - cos(theta_2)*sin(i)*del_RAAN);

rho_2H_new = vpa([x2,y2,z2])/1000


% # adjust the return matrix values as needed
% def result():
%     rho_0_H = [-1.3393, 2.5279, 6.4976] # km
%     rho_t_H = [-2.2602, 7.0763, 10.1773] # km
% 
%     return rho_0_H, rho_t_H
% result()

%% Quiz 17 - Relative Motion Trajectories w/ J2 Perturbations:
% See page 70 on my notes!
J2 = 1.08263e-3; r_eq = 6378e3; 
a = 8000e3; e = 0.1; i = 80*pi/180;
del_a = 10e3; del_e = 0.05; del_i = 1*pi/180;

eta = sqrt(1 - e^2);
n = sqrt(mu/a^3);
epsilon = 3*J2*(r_eq/(a*(1-e^2)))^2;
del_kappa_RAAN = 7/4*cos(i)*del_a/a - 2*e/eta^2 * cos(i)*del_e + 1/2*sin(i)*del_i;

RAAN_Dot = epsilon * n * del_kappa_RAAN
t_period = 2*pi*sqrt(a^3/mu);

del_RAAN = RAAN_Dot * t_period
