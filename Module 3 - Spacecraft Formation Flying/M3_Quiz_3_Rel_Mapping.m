%% 
clear all; close all; clc;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Relative Motion Mapping (Orbit-Hill <-> Inertial)
% Chief-Deputy mapping:
% r_d,N = r_c,N + [NH] * rho,H

%% Q1 - Map the inertial chief and deputy velocity vectors to the equivalent Hill frame relative position vector and velocity.
global mu;
mu = 3.986e14; % [km^3/s^2]

r_cN = 1e3 * [-4893.268,3864.478,3169.646]';
v_cN = 1e3 * [-3.91386,-6.257673,1.59797]';

r_dN = 1e3 * [-4892.98,3863.073,3170.619]';
v_dN = 1e3 * [-3.913302,-6.258661,1.598199]';

[rho_H,rhoP_H] = rv2hill(r_cN, v_cN, r_dN, v_dN);
vpa(rho_H)/1e3
vpa(rhoP_H)/1e3

% # adjust the return matrix values as needed
% def result():
%     rho_H = [-0.53680960588475545591791160404682,1.2211950820166237008379539474845,1.1064429577200378389534307643771]  # km
%     rhoP_H = [0.00048602514641049554544594002436497,0.0011580569652054402141061473230366,0.00055886342117941112395129721335252] # km/s
%     return rho_H, rhoP_H
% 
% result()

%% Q2 - Reverse mapping:
rho_H2 = 1e3*[-0.537,1.221,1.106]';
rhoP_H2 = 1e3*[0.000486,0.001158,0.0005590]';
[r_dN2, v_dN2] = hill2rv(r_cN, v_cN, rho_H2, rhoP_H2);

vpa(r_dN2)/1e3
vpa(v_dN2)/1e3

% # adjust the return matrix values as needed
% def result():
%     rd_N = [-4892.9799839025214314460754394531,3863.0730949881165288388729095459,3170.6184888869174756109714508057]  # km
%     vd_N = [-3.913301926561904110712930560112,-6.2586606917541303118923678994179,1.5981991468557403095474001020193] # km/s
%     return rd_N, vd_N
% 
% result()

%% Q3 - Rectilinear -> Curvilinear: (SEE iPad for derivation)
x = 10; y = 500; r = 7e3; xDot = 0.1; yDot = -0.1; rDot = 0.05;
rd = sqrt((r+x)^2 + y^2);
dr = rd - r
% dtheta = acos((r+x)/rd)
dtheta = atan2(y,(r+x));
s = r*dtheta

rdDot = ((r+x)*(rDot + xDot) + y*yDot) / rd;
drDot = rdDot - rDot
dthetaDot = ( (r+x)*yDot - (rDot + xDot)*y ) / rd^2;
sDot = rDot*dtheta + r*dthetaDot

% # adjust the return matrix values as needed
% def result():
%     delta_r = [27.8090]  # km
%     delta_r_dot = [0.0925] # km/s
%     s = [498.4426] # km
%     s_dot = [-0.1064] # km/s
%     return delta_r, delta_r_dot, s, s_dot
% 
% result()

%% Q4 - Rectilinear -> Curvilinear:
dr2 = 10; s2 = 500; r2 = 7e3; drDot2 = 0.1; sDot2 = -0.1; rDot2 = 0.05;
rd2 = r2 + dr2; 
rdDot2 = rDot2 + drDot2;
dtheta2 = s2/r2;
dthetaDot2 = (sDot2 - rDot2*dtheta2)/r2;

x2 = rd2*cos(dtheta2) - r2
y2 = (r2 + x2)*tan(dtheta2)

A = [r2 + x2, y2; -y2, r2+x2];
B = [rd2*rdDot2 - rDot2*(r2+x2); dthetaDot2*rd2^2 + rDot2*y2];
sol = A\B

% # adjust the return matrix values as needed
% def result():
%     x = [-7.8751]  # km
%     x_dot = [0.1070] # km/s
%     y = [500.2886] # km
%     y_dot = [-0.0927] # km/s
%     return x, x_dot, y, y_dot
% 
% result()


%% Function
function [rho_H, rhoP_H] = rv2hill(r_cN, v_cN, r_dN, v_dN)
    % global mu;
    [ON,~,omega_HN_H] = N2H(r_cN, v_cN);

    rho_H = ON * (r_dN - r_cN);
    rhoP_H = ON * (v_dN - v_cN) - cross(omega_HN_H, rho_H);
end

function [r_dN, v_dN] = hill2rv(r_cN, v_cN, rho_H, rhoP_H)
    [ON,~,omega_HN_H] = N2H(r_cN, v_cN);
    r_dN = r_cN + ON' * rho_H;
    v_dN = v_cN + ON' * (rhoP_H + cross(omega_HN_H,rho_H));
end

function [ON,h,omega_HN_H] = N2H(r_cN, v_cN)
    h_vec = cross(r_cN, v_cN);
    h = norm(h_vec);

    o_cap_r = r_cN/norm(r_cN);
    o_cap_h = h_vec/h;
    o_cap_theta = cross(o_cap_h, o_cap_r);

    ON = [o_cap_r';o_cap_theta';o_cap_h'];
    fDot = h/(norm(r_cN)^2);
    omega_HN_H = [0,0,fDot]';
end