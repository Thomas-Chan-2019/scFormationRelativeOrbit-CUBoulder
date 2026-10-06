%% 
clear all; close all; clc;
import astro.*   % COE2RV / RV2COE from lib/+astro (run `startup` if not found)
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Quiz 4 - Linearizing
% Q2
A = 1/16 * [-10, 5*sqrt(2)*(6-5*sqrt(3)), 5*sqrt(2)*(6+5*sqrt(3)); ...
            5*sqrt(2)*(6-5*sqrt(3)), 5-20*sqrt(3), 115; ...
            5*sqrt(2)*(6+5*sqrt(3)), 115, 5+20*sqrt(3)]

[V,D,W] = eig(A)

ans = V(:,1)