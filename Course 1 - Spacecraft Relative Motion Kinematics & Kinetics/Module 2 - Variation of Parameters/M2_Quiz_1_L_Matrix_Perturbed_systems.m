%% 
clear all; close all; clc;

%% Q4: Lagrange Matrix:
tf = 10;
e1_0 = 1;
e2_0 = 0;
x0 = e1_0; 
xDot_0 = e2_0;
s0 = [x0, xDot_0, e1_0, e2_0]';

[t, s] = ode45(@odefunc,[0 10],s0);

%%
figure;
plot(t, s);
legend('$x$','$\dot{x}$','$e_1$','$e_2$');
xlabel('$t$ [s]');


%% functions
function ds = odefunc(t,s)
% Dynamical systems: x'' + 3x + x^3 = 0 => x'' = -3x - x^3 = -(sqrt(3))^2*x + (-x^3) 
% -> Using Perturbed Spring Mass form:
% x'' = -w^2 * x + a_d(t,x,x',...), w = frequency
% So: x'' = -(sqrt(3))^2*x + (-x^3), where: a_d = -x^3

% Solution for Perturbed Spring Mass systems:
% Assuming e1 = x0, e2 = xDot_0,
% f(t, e(t)) := x(t) = e1(t) * cos(w*t) + e2(t)/w * sin(w*t);
% xDot(t) = -w*e1(t) * cos(w*t) + e2(t)*sin(w*t);

% Solving for [L]: [L] = [df/de1 df/de2; d^2f/dtde1 d^2f/dtde2] 
%                      = [cos(w*t), 1/w * sin(w*t); -w*sin(w*t), cos(w*t)] 
% Therefore: [L]*eDot = [L]*[e1Dot; e2_Dot] = [0; a_d]
%            e1Dot = dx0/dt = - 1/w * sin(w*t) * a_d
%            e2Dot = dxDot_0/dt = cos(w*t) * a_d

% To solve the system above, arrange state s: 
% s = [x, xDot, e1, e2]'
    x = s(1);
    xDot = s(2);
    e1 = s(3);
    e2 = s(4);

    xDDot = -3*x -x^3;
    e1Dot = 1/sqrt(3) * x^3 * sin(sqrt(3)*t);
    e2Dot = - x^3 * cos(sqrt(3)*t);
    
    ds = [xDot, xDDot, e1Dot, e2Dot]';
end 