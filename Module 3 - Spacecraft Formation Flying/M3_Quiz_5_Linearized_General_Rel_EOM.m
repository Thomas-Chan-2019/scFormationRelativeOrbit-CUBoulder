%% 
clear all; close all; clc;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Q2 - Deputy Large Rel.Motion w/o perturbations - full nonlinear simulations
global mu; mu = 3.986e14;

tf = 2000;
r_cN0 = 1e3 * [-6685.20926, 601.51244, 3346.06634]';
v_cN0 = 1e3 * [-1.74294,-6.70242,-2.27739]';

rho_H0 = 1e3 * [-81.22301, 248.14201, 94.95904]';
rho_P_H0 = 1e3 * [0.47884, 0.14857, 0.13577]';

%% Task 1 - Map deputy Hill to inertial state for Inetial integration:
% [r_dN0, v_dN0] = hill2rv(r_cN0,v_cN0,rho_H0,rho_P_H0);
% 
% s0_deputy_N = [r_cN0; v_cN0; r_dN0; v_dN0];
% 
% [t, s_deputy_N] = ode45(@InertialODE, [0 tf], s0_deputy_N);
% 
% figure; 
% subplot(2,2,1); plot(t,s_deputy_N(:,7:9)); title('Inertial Integration - $r_{d/N,1}$');
% subplot(2,2,2); plot(t,s_deputy_N(:,10:end)); title('Inertial Integration - $v_{d/N,1}$');
% 
% %% 
% r_dN_f = s_deputy_N(end,7:9)';
% v_dN_f = s_deputy_N(end,10:end)';
% 
% vpa(r_dN_f/1e3)
% vpa(v_dN_f/1e3)
% 
% % Get chief states:
% % [t, s_chief_N] = ode45(@InertialODE, [0 tf], [r_cN0;v_cN0]);
% r_cN_f = s_deputy_N(end,1:3)';
% v_cN_f = s_deputy_N(end,4:6)';
% 
% [rho_H_f,rho_P_H_f] = rv2hill(r_cN_f,v_cN_f,r_dN_f,v_dN_f)

%% Task 2 - Use nonlinear rel. motion equation for integration:
s0_deputy_Rel = [r_cN0; v_cN0; rho_H0; rho_P_H0];

[t2, s_deputy_Rel] = ode45(@linRelMotionODE, [0 tf], s0_deputy_Rel);
s_deputy_Rel = s_deputy_Rel'; % transposing for correct dimension.

[r_dN_f_rel, v_dN_f_rel] = hill2rv(s_deputy_Rel(1:3,end), s_deputy_Rel(4:6,end), s_deputy_Rel(7:9,end), s_deputy_Rel(10:end,end));
vpa(r_dN_f_rel/1e3)
vpa(v_dN_f_rel/1e3)
s_deputy_Rel(7:9,end)/1e3
s_deputy_Rel(10:end,end)/1e3

r_dN_rel = [];
v_dN_rel = [];
for i = 1:length(t2)
    [r_dN_rel_temp, v_dN_rel_temp] = hill2rv(s_deputy_Rel(1:3,i), s_deputy_Rel(4:6,i), s_deputy_Rel(7:9,i), s_deputy_Rel(10:end,i));
    r_dN_rel(:,i) = r_dN_rel_temp;
    v_dN_rel(:,i) = v_dN_rel_temp;
end

figure; hold on;
subplot(2,2,3); plot(t2,r_dN_rel); title('Relative Integration - $r_{d,N,2}$');
subplot(2,2,4); plot(t2,v_dN_rel); title('Relative Integration - $v_{d,N,2}$');

%% Check r_cN:
% figure; hold on;
% subplot(2,1,1); plot(t2,s_deputy_Rel(1:3,:)); title('Relative Integration - $r_{d,N,2}$');
% subplot(2,1,2); plot(t2,s_deputy_Rel(4:6,:)); title('Relative Integration - $v_{d,N,2}$');

%% Check Deputy Final Inertial state:
% r_dN_f - r_dN_f_rel
% v_dN_f - v_dN_f_rel

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

function ds = linRelMotionODE(t,s)
% Given  
    global mu;
    r_cN_vec = s(1:3); 
    r_cN = norm(r_cN_vec);

    r_cNDot_vec = s(4:6); 
    % r_cNDot = norm(r_cNDot_vec); % Wrong Frame! We are at Hill Orbit
    % frame so cannot directly take norm.
    r_cNDot = dot(r_cNDot_vec, r_cN_vec/r_cN); % Correct! VERY IMPORTANT!
    
    rho_H = s(7:9); 
    x = rho_H(1); y = rho_H(2); z = rho_H(3);

    rho_P_H = s(10:end);
    xDot = rho_P_H(1); yDot = rho_P_H(2); zDot = rho_P_H(3);
    
    % [ON,~,~] = N2H(r_cN_vec, r_cNDot_vec);
    % r_dN_vec = r_cN_vec + ON' * rho_H;
    % r_dN_norm = norm(r_dN_vec)

    r_dN_norm = norm([r_cN + x, y, z]);

    % h = r_cN^2 * fDot; f always refer to chief!
    h = norm(cross(r_cN_vec, r_cNDot_vec));
    fDot = h/r_cN^2;
    p = h^2/mu;

    r_cNDDot_vec = -mu/r_cN^3 * r_cN_vec;
    rho_PP_H = [ x*fDot^2*(1 + 2*r_cN/p) + 2*fDot*(yDot - y*r_cNDot/r_cN); ...
               - 2*fDot*(xDot - x*r_cNDot/r_cN) + y*fDot^2*(1 - r_cN/p); ...
               - r_cN/p * fDot^2 * z];

    ds = [r_cNDot_vec; r_cNDDot_vec; rho_P_H; rho_PP_H];
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