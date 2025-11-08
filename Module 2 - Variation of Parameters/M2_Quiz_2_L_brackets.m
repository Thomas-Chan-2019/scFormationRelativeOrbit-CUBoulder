%% 
clear all; close all; clc;

%% Q2: Lagrange Brackets:
tf = 10;
e1_0 = 1;
e2_0 = 0;
y0 = e1_0; 
yDot_0 = e2_0;
s0 = [y0, yDot_0, e1_0, e2_0]';

[t, s] = ode45(@odefunc,[0 tf],s0);

%%
figure;
plot(t, s);
legend('$x$','$\dot{x}$','$e_1$','$e_2$');
xlabel('$t$ [s]');


%% functions
function ds = odefunc(t,s)
% Given potential function and consider only vertical motion y(t): 
% V(y) = mgy + k/2 * y^2, where m = 1, g = 9.81 m/s^2, k = 3.
%      = 9.81*y + 1.5*y^2 
% Dis. potential: R(y) := -g*y
% EOM: m*y'' = -dV/dy = -9.81 - 3*y 
% Again using Perturbed Spring Mass systems:
% y'' = -w^2*y + a_d(t,y,y',...)
%     = -(sqrt(3))^2*y + (a_d == -9.81)
% => See GoodNote!

% e1Dot = dx0/dt = - 1/w * sin(w*t) * a_d
% e2Dot = dxDot_0/dt = cos(w*t) * a_d
    y = s(1);
    yDot = s(2);
    e1 = s(3);
    e2 = s(4);

    g = 9.81;
    a_d = -g;

    yDDot = -3*y + a_d;
    e1Dot = g/sqrt(3) * sin(sqrt(3)*t);
    e2Dot = -g * cos(sqrt(3)*t);
    
    ds = [yDot, yDDot, e1Dot, e2Dot]';
end 