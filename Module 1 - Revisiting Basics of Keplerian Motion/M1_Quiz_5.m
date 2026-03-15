close all;
clear all;
%% Q3
r0 = [2466.69,5941.54,3282.71]*1e3;
v0 = [-6.80822,1.04998,3.61939]*1e3;

s0 = [r0,v0]';

h = norm(cross(r0,v0));

global mu
mu = 3.986e14; % mu Earth

[t, s] = ode45(@q3ODE,[0 3600],s0);

rf = s(end,1:3);
vf = s(end,4:6);
assert(abs(h - norm(cross(rf,vf)))/h <= 0.001)
% vpa(norm(cross(rf,vf)))
% vpa(h)

% Python result submitted:
% # adjust the return matrix values as needed
% def result():
%     rVec_N = [-3436.22, -7152.73, -3756.10]  # km
%     vVec_N = [5.57478, -0.92151, -3.00854]  # km/sec
% 
%     return rVec_N, vVec_N
% 
% result()

%% Q4
r1_0 = [-6685.20926, 601.51244, 3346.06634]*1e3;
v1_0 = [-1.74294, -6.70242, -2.27739]*1e3;
r2_0 = [-6685.21657, 592.52839, 3345.6716]*1e3;
v2_0 = [-1.74283, -6.70475, -2.27334]*1e3;

s0 = [r1_0, v1_0, r2_0, v2_0]';
[t, s] = ode45(@q4ODE,[0 10104848],s0)

s1_f = s(end,1:6) / 1e3;
s2_f = s(end,7:end) / 1e3;

r1_f = vpa(s1_f(1:3))
v1_f = vpa(s1_f(4:end))
r2_f = vpa(s2_f(1:3))
v2_f = vpa(s2_f(4:end))

% Python result submitted:
% # adjust the return matrix values as needed
% def result():
%     r1 = [1737.42736, 6888.84178, 2381.18927] # km
%     v1 = [-6.52127, 0.53217, 3.22993] # km/sec
%     r2 = [1745.29067, 6880.45424, 2368.44662] # km
%     v2 = [-6.52590, 0.54282, 3.24066] # km/sec
% 
%     return r1, v1, r2, v2
% 
% result()

figure; 
% subplot(1,2,1); 
hold on;
plot3(s(:,1),s(:,2),s(:,3));
plot3(0,0,0,'rx');
% 
% subplot(1,2,2); hold on;
plot3(s(:,7),s(:,8),s(:,9));
plot3(0,0,0,'rx');
legend;


%% Funcions:
function ds = q3ODE(t,s)
    global mu;
    r = s(1:3);
    r_Dot = s(4:6);

    r_Ddot = (-mu/norm(r)^3) * r;

    ds = [r_Dot; r_Ddot];
end 

function ds = q4ODE(t,s)
    global mu;
    s1 = s(1:6);
    s2 = s(7:end);
    
    r1 = s1(1:3);
    r1_Dot = s1(4:6);

    r2 = s2(1:3);
    r2_Dot = s2(4:6);

    r1_Ddot = (-mu/norm(r1)^3) * r1 + perturbJ2(r1);
    r2_Ddot = (-mu/norm(r2)^3) * r2 + perturbJ2(r2);

    ds = [r1_Dot; r1_Ddot; r2_Dot; r2_Ddot];
end 

function a_J2 = perturbJ2(r)
    global mu;
    x = r(1); y = r(2); z = r(3);
    
    r_norm = norm(r);
    J2 = 1082.63e-6;
    r_eq = 6378.14e3;
    
    a_vec = [(1 - 5*(z/r_norm)^2) * x/r_norm; (1 - 5*(z/r_norm)^2) * y/r_norm;(3 - 5*(z/r_norm)^2) * z/r_norm];
    
    a_J2 = -3/2 * J2 * (mu/r_norm^2) * (r_eq/r_norm)^2 * a_vec;
end