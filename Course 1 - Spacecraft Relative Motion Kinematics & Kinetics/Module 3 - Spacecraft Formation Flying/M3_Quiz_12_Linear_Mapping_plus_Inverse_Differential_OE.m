%% 
clear all; close all; clc;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Q3 - Linear Mapping: Cartesian => Differential OE:
global mu; 
mu = 3.986e5; % [km3/s2], Use KM !

OED = zeros(6,1);
OE_c = [7500, deg2rad(13), deg2rad(22), 0.00707107, 0.00707107, deg2rad(70)]';
% OE_c = [7500e3, deg2rad(13), deg2rad(22), 0.00707107, 0.00707107, deg2rad(70)]';
[~, A] = cart2OED(OE_c, OED);
A
vpa(A)

% ANS in Python:
% # adjust the return matrix values as needed
% def result():
%     row1 = [0.99149167109588420832011479433277, -39.074432374736659312475239858031, 0, -7289.8431830407562301843427121639, -1763.8856724042250334605341777205, 0]  #
%     row2 = [0, 7436.1875332191311827045865356922, 0, 0, 0, 6894.7130203316046390682458877563]  #
%     row3 = [0, 0, 1672.7782258963482036051573231816, 0, 0, -2714.2489820199371024500578641891]  #
%     row4 = [0.0000025755993950768936348030488381333, 0.061827309673782818166820618444035, 0, 1.6397415779907305033447073583375, -7.1039577028510558776019934157375, 0]  #
%     row5 = [-0.0014704735768308093823292237445344, 0.038633990926153403966392119173179, 0, 14.363351905977138400771764281671, 3.4360125001578176728855851251865, -0.035820812624209205643488473924663]  #
%     row6 = [0, 0, 7.155236414615712803310998424422, 0, 0, 0.63367204340796068340324609380332]  #
% 
%     return row1, row2, row3, row4, row5, row6
% 
% result()

%% Q4 - Inverse Mapping: Differential OE => Cartesian:
dX = zeros(6,1);
[~, invA] = OED2cart(OE_c, dX);
invA
vpa(invA)

% ANS in Python:
% # adjust the return matrix values as needed
% def result():
%     row1 = [4.0861986038467668080897965410259, 0.010781065812809175266662187198108, 0, -10.903973856478321735608005838003, 2075.1163747298814996611326932907, 0]  #
%     row2 = [0, 0.0001344775122376586009177029668038, 0.00032391933196776500399755027714832, 0, 0, -0.075727086299449589135157623331907]  #
%     row3 = [0, 0, 0.000030939377185482049318829278483634, 0, 0, 0.13252466146113484057877940358594]  #
%     row4 = [0.00039610559181547905878081894393006, 0.00000043697374546006309379650596989852, -0.0000022904562706973041036860594399371, 0.030855195447547667519927472312702, 0.26713493664325793641722839311115, 0.0005354715281194489260099822480754]  #
%     row5 = [0.000092910916201232939893098206063371, 0.0000012751588477432694964985903166022, 0.0000022904562706973041036860594399371, -0.13364853467737503045498215215048, 0.062412664823355591681153242689106, -0.0005354715281194489260099822480754]  #
%     row6 = [0, 0, -0.00034935825335214820523688228348647, 0, 0, 0.081674293535700007162780877933983]  #
% 
%     return row1, row2, row3, row4, row5, row6
% 
% result()

%% functions
function [dX,A] = cart2OED(OE_c, OED) % CAUTION -> THIS IS NON-SI UNIT -> USE KM !
    global mu;

    a = OE_c(1);
    theta = OE_c(2);
    i = OE_c(3);
    q1 = OE_c(4);
    q2 = OE_c(5);
    RAAN = OE_c(6);

    e = sqrt(q1^2 + q2^2); 
    omega = acos(q1/e);
    f = theta - omega;
    if f<0
        f = f + 2*pi;
    end
    p = a * (1 - q1^2 - q2^2);
    h = sqrt(p * mu);

    R = a * (1 - q1^2 - q2^2) / (1 + e*cos(f));
    Vr = h/p * (q1*sin(theta) - q2*cos(theta));
    Vt = h/p * (1 + q1*cos(theta) + q2*sin(theta));

    A = [R/a, Vr/Vt*R,  0, -R/p * (2*a*q1 + R*cos(theta)), -R/p * (2*a*q2 + R*sin(theta)),  0;...
          0,  R,  0,  0,  0, R*cos(i);...
          0,  0, R*sin(theta),  0,  0, -R*cos(theta)*sin(i);...
         -Vr/(2*a), (1/R - 1/p)*h,  0, (Vr*a*q1 + h*sin(theta))/p, (Vr*a*q2 - h*cos(theta))/p,  0;...
         -3*Vt/(2*a), -Vr,  0, (3*Vt*a*q1 + 2*h*cos(theta))/p, (3*Vt*a*q2 + 2*h*sin(theta))/p, Vr*cos(i);...
          0,  0, Vt*cos(theta) + Vr*sin(theta),  0,  0, (Vt*sin(theta) - Vr*cos(theta))*sin(i)
        ];
    
    dX = A * OED;
end 

function [OED, invA] = OED2cart(OE_c, dX)
    global mu;
    a = OE_c(1);
    theta = OE_c(2);
    i = OE_c(3);
    q1 = OE_c(4);
    q2 = OE_c(5);
    RAAN = OE_c(6);

    e = sqrt(q1^2 + q2^2); 
    omega = acos(q1/e);
    f = theta - omega;
    if f<0
        f = f + 2*pi;
    end
    p = a * (1 - q1^2 - q2^2);
    h = sqrt(p * mu);

    R = a * (1 - q1^2 - q2^2) / (1 + e*cos(f));
    Vr = h/p * (q1*sin(theta) - q2*cos(theta));
    Vt = h/p * (1 + q1*cos(theta) + q2*sin(theta));

    % Extra intermediate params:
    alpha = a/R;
    nu = Vr/Vt;
    rho = R/p;
    kappa1 = alpha * (1/rho - 1);
    kappa2 = alpha * nu^2 / rho;

    invA = [2*alpha*(2 + 3*kappa1 + 2*kappa2), -2*alpha*nu*(1 + 2*kappa1 + kappa2),  0, 2*alpha^2*nu*p/Vt, 2*a*(1 + 2*kappa1 + kappa2)/Vt,  0;...
             0, 1/R, cot(i)*(cos(theta) + nu*sin(theta))/R,  0,  0, -sin(theta)*cot(i)/Vt;...
             0,  0, (sin(theta) - nu*cos(theta))/R,  0,  0, cos(theta)/Vt;...
            (3*cos(theta) + 2*nu*sin(theta))/(rho*R), -(nu^2*sin(theta)/rho + q1*sin(2*theta) - q2*cos(2*theta))/R, -q2*cot(i)*(cos(theta) + nu*sin(theta))/R,  sin(theta)/(rho*Vt), (2*cos(theta) + nu*sin(theta))/(rho*Vt),  q2*cot(i)*sin(theta)/Vt;...
            (3*sin(theta) - 2*nu*cos(theta))/(rho*R),  (nu^2*cos(theta)/rho + q2*sin(2*theta) + q1*cos(2*theta))/R,  q1*cot(i)*(cos(theta) + nu*sin(theta))/R, -cos(theta)/(rho*Vt), (2*sin(theta) - nu*cos(theta))/(rho*Vt), -q1*cot(i)*sin(theta)/Vt;...
             0,  0, -(cos(theta) + nu*sin(theta))/(R*sin(i)),  0,  0, sin(theta)/(Vt*sin(i))
           ];

    OED = invA * dX;

end