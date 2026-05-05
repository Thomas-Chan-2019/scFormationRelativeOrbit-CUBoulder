%% 
clear all; close all; clc;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Q1 - Integrate for CW
global mu a_c n; 
mu = 3.986e14; 
% a_c = 6800e3; n = sqrt(mu/a_c^3);
% 
% tf = 1300;
% 
% rho_H0 = 1e3 * [1.299038, -1.0000, 0.3213938]';
% rho_P_H0 = 1e3 * [-0.000844437, -0.00292521, -0.000431250]';
% 
% % Find rho, rho':
% s0_deputy_Rel = [rho_H0; rho_P_H0];
% 
% [~, s_deputy_Rel] = ode45(@CWH_ODE, [0 tf], s0_deputy_Rel);
% s_deputy_Rel = s_deputy_Rel'; % transposing for correct dimension.
% 
% disp("Q1 - [km] rho_H = ");
% disp((s_deputy_Rel(1:3,end)/1e3)');
% disp("Q1 - [km/s] rho_PH = ");
% disp((s_deputy_Rel(4:end,end)/1e3)');

%% Q2: Compare CWH with Inertial Full-Nonlinear solution:
% Leader follower scenario w/ y_off = 200 km along track
a_c = 7500e3; 
n = sqrt(mu/a_c^3);
rho_H0 = 1e3 * [0, 200, 0]';
rho_P_H0 = 1e3 * [0, 0, 0]';

% theta = 200e3 / a_c; 
% r_cN0 = [a_c; 0; 0]; 
% v_cN0 = [0; sqrt(mu/a_c); 0]; 
% r_dN0 = a_c * [cos(theta); sin(theta); 0]; 
% v_dN0 = sqrt(mu/a_c) * [-sin(theta); cos(theta); 0]; 
% [rho_H0, rho_P_H0] = rv2hill(r_cN0, v_cN0, r_dN0, v_dN0);

tf = 2000;

s0_deputy_Rel = [rho_H0; rho_P_H0];

[t2, s_deputy_Rel] = ode45(@CWH_ODE, [0 tf], s0_deputy_Rel);
s_deputy_Rel = s_deputy_Rel'; % transposing for correct dimension.

rho_H_f = (s_deputy_Rel(1:3,end));
rho_PH_f = (s_deputy_Rel(4:end,end));

disp("Q2 - [km] rho_H = ");
disp(rho_H_f/1e3);
disp("Q2 - [km/s] rho_PH = ");
disp(rho_PH_f/1e3);

% Inertial integrations:
% r_cN0 = 1e3 * [-6685.20926, 601.51244, 3346.06634]';
% v_cN0 = 1e3 * [-1.74294,-6.70242,-2.27739]';
r_cN0 = [a_c, 0, 0]';
v_cN0 = [0, sqrt(mu/a_c), 0]';

[r_dN0, v_dN0] = hill2rv(r_cN0,v_cN0,rho_H0,rho_P_H0);

s0_deputy_N = [r_cN0; v_cN0; r_dN0; v_dN0];
[t, s_deputy_N] = ode45(@InertialODE, [0 tf], s0_deputy_N);

r_cN_f = s_deputy_N(end,1:3)';
v_cN_f = s_deputy_N(end,4:6)';
r_dN_f = s_deputy_N(end,7:9)';
v_dN_f = s_deputy_N(end,10:end)';

[rho_H_f_N,rho_PH_f_N] = rv2hill(r_cN_f,v_cN_f,r_dN_f,v_dN_f);


disp("Q2 - Inertial - [km] rho_H = ");
disp(vpa(rho_H_f_N/1e3));
% disp("Q2 - Inertial - [km/s] rho_PH = ");
% disp(rho_PH_f_N/1e3);

% Asked: TOTAL SEPARATION DISPLACEMENT (WITH UNITS FUKKKK) !!!
err_pos = norm(rho_H_f_N) - norm(rho_H0); 
disp("Error distance [km]: ");
disp(vpa(err_pos)/1e3);
% disp("Error velocity between Nonlinear & CWH:");
% disp((rho_PH_f_N - rho_PH_f)/1e3);

% Nonlinear relative motion method:
s0_deputy_Rel = [r_cN0; v_cN0; rho_H0; rho_P_H0];

[t2, s_deputy_Rel] = ode45(@relMotionODE, [0 tf], s0_deputy_Rel);
s_deputy_Rel = s_deputy_Rel'; % transposing for correct dimension.

[r_dN_f_rel, v_dN_f_rel] = hill2rv(s_deputy_Rel(1:3,end), s_deputy_Rel(4:6,end), s_deputy_Rel(7:9,end), s_deputy_Rel(10:end,end));

rho_H_f_rel = s_deputy_Rel(7:9,end);
disp("Q2 - Rel - [km] rho_H = ");
disp(vpa(rho_H_f_rel)/1e3);

%% Function
function ds = InertialODE(t,s)
% Given Keplerian Unperturbed motion: 
% r'' = - mu/|r|^3 * r
    global mu;
    r_cN = s(1:3);
    r_cNDot = s(4:6);
    
    r_dN = s(7:9);
    r_dNDot = s(10:end);
    
    r_cNDDot = ( -mu/norm(r_cN)^3 ) * r_cN;
    r_dNDDot = ( -mu/norm(r_dN)^3 ) * r_dN;
    

    ds = [r_cNDot; r_cNDDot; r_dNDot; r_dNDDot];
end 

function ds = relMotionODE(t,s)
% Given  
    global mu;
    r_cN_vec = s(1:3); 
    r_cN = norm(r_cN_vec);

    r_cNDot_vec = s(4:6); 
    r_cNDot = dot(r_cNDot_vec, r_cN_vec/r_cN);
    
    rho_H = s(7:9); 
    x = rho_H(1); y = rho_H(2); z = rho_H(3);

    rho_P_H = s(10:end);
    xDot = rho_P_H(1); yDot = rho_P_H(2); zDot = rho_P_H(3);
 
    % r_dN_norm = norm([r_cN + x, y, z]);
    [ON,~,~] = N2H(r_cN_vec, r_cNDot_vec);
    r_dN_vec = r_cN_vec + ON' * rho_H;
    r_dN_norm = norm(r_dN_vec);

    % h = r_cN^2 * fDot; f always refer to chief!
    h = norm(cross(r_cN_vec, r_cNDot_vec));
    fDot = h/r_cN^2;

    r_cNDDot_vec = -mu/r_cN^3 * r_cN_vec;
    rho_PP_H = [-mu/r_dN_norm^3 * (r_cN + x) + 2*fDot*(yDot - y*r_cNDot/r_cN) + x*fDot^2 + mu/r_cN^2; ...
                -mu/r_dN_norm^3 * y - 2*fDot*(xDot - x*r_cNDot/r_cN) + y*fDot^2; ...
                -mu/r_dN_norm^3 * z];

    ds = [r_cNDot_vec; r_cNDDot_vec; rho_P_H; rho_PP_H];
end

function ds = CWH_ODE(t,s)
% Given  
    global n;
    
    rho_H = s(1:3); 
    x = rho_H(1); y = rho_H(2); z = rho_H(3);

    rho_P_H = s(4:end);
    xDot = rho_P_H(1); yDot = rho_P_H(2); zDot = rho_P_H(3);
    
    rho_PP_H = [  2*n*yDot + 3*n^2*x; ...
                - 2*n*xDot; ...
                - n^2*z];

    ds = [rho_P_H; rho_PP_H];
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