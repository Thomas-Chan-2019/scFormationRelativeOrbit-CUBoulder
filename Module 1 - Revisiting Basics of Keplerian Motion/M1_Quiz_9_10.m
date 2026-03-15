clear all; close all;
%% Quiz 9 - Fundamental Integrals & Kepler's Equation (i.e., Time of Flight)
a = 7500e3; 
e = 0.05;
f0 = 25*pi/180; % [rad], initial true anomaly
tf = 3600; % 1 hour 

% % Check answer from Arizona State University - Spacecraft Dynamics and Control course.
% a = 25512e3; 
% e = 0.625;
% f0 = 0*pi/180; % [rad], initial true anomaly
% tf = 4*3600; % 1 hour 

mu = 3.986e14; % Earth orbit
n = sqrt(mu/a^3); % [s^-1], Mean motion

% Quiz 9 Steps:
E0 = f2E(f0,e)
M0 = E2M(E0,e)

Mf = M0 + n*(tf - 0)
Ef = M2E(Mf,e)

% EccentricAnomaly(Mf,e)

ff = E2f(Ef,e)

%% Quiz 10 - Coordinate Transformation:
% Q6 - COE to RV:
a = 8000e3; 
e = 0.1;
i = 30*pi/180;
RAAN = 145*pi/180;
omega = 120*pi/180;
M0 = 10*pi/180;

tf = 3600; % 1 hour later
mu = 3.986e14; % Earth orbit
n = sqrt(mu/a^3); % [s^-1], Mean motion

% Quiz 10 Steps:
Mf = M0 + n*(tf - 0);
Ef = M2E(Mf,e);
ff = E2f(Ef,e);

theta = omega + ff;
r_peri = a*(1-e^2)/(1 + e*cos(ff));
% x_peri = r_peri*cos(ff);
% y_peri = r_peri*sin(ff);

h = sqrt(mu*a*(1-e^2));

r_vec = r_peri * [cos(RAAN)*cos(theta) - sin(RAAN)*sin(theta)*cos(i);...
           sin(RAAN)*cos(theta) + cos(RAAN)*sin(theta)*cos(i);...
           sin(theta)*sin(i)]

v_vec = -mu/h * [cos(RAAN)*(sin(theta) + e*sin(omega)) + sin(RAAN)*(cos(theta) + e*cos(omega))*cos(i);...
                 sin(RAAN)*(sin(theta) + e*sin(omega)) - cos(RAAN)*(cos(theta) + e*cos(omega))*cos(i);...
                 -(cos(theta) + e*cos(omega))*sin(i)]

% Check
% assert((norm(cross(r_vec,v_vec)) - h)/h <= 1e-6)

% Q7 - RV to COE:
r_vec = [-820.865, -1905.95, -7445.9]' * 1e3;
v_vec = [-6.75764, -1.85916, 0.930651]' * 1e3;

[a,e,i,RAAN,omega,f] = RV2COE(r_vec,v_vec,mu)

%% Functions
function M = E2M(E,e)
    M = E - e*sin(E);
end

function E = f2E(f,e) % f in rad.
    E = 2 * atan(sqrt((1-e)/(1+e)) * tan(f/2));
end

function f = E2f(E,e) % E in rad.
    f = 2 * atan(sqrt((1+e)/(1-e)) * tan(E/2));
end

function E = M2E(M,e) % M in rad.
% Take reference if needed, from function `EccentricAnomaly(M_k,e)` from KTH GNSS course path: 
% C:\Users\chany\OneDrive\文件\KTH-Academics-New\Sem 2 Year 1 (Spring 2023)\AH2923 Global Navigation Satellite Systems\ASM 3&4\asm_4_receiver_pos.m
    E0 = M;
    Ek = E0;
    
    Ef = 0; % Final E.
    while true
      % Ek_1 = Ek - (M - Ek + e*sin(Ek))/(e*cos(Ek) - 1);
      Ek_1 = M + e*sin(Ek);
      
      % Break condition as E solution converges within 10^-13:
      if abs(Ek_1 - Ek) <= 1e-13
        Ef = Ek_1;
        break;
      end
      Ek = Ek_1;
    end

    E = Ef;
end

function [a,e,i,RAAN,omega,f] = RV2COE(r,v,mu)
    % General reused terms:
    hVec = cross(r,v);
    h = norm(hVec);
    
    iCap_r = r/norm(r);

    % Semimajor axis a:
    a_inv = 2/norm(r) - norm(v)^2/mu;
    if (a_inv ~= 0)
        a = 1/a_inv; 
    else
        error("Semimajor axis a goes to infinity. Check the orbit conic section!");
    end

    % Eccentricity e:
    eVec = cross(v,hVec)/mu - iCap_r;
    e = norm(eVec);

    % Construct PN (perifocal to inertial frame DCM):
    iCap_e = eVec/e; iCap_h = hVec/h; iCap_theta = cross(iCap_h,iCap_e);
    PN = [iCap_e,iCap_theta,iCap_h]';
    
    RAAN = atan2(PN(3,1),-PN(3,2));
    i = acos(PN(3,3));
    omega = atan2(PN(1,3),PN(2,3));

    f = atan2(dot(cross(iCap_e,iCap_r),iCap_h), dot(iCap_e,iCap_r));
end