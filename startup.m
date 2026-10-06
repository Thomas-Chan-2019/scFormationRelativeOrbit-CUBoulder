% startup.m  Register this repo's shared MATLAB libraries on the path.
%
% MATLAB runs a startup.m automatically at launch IF it is on the path at
% startup (i.e. when MATLAB's startup folder is this repo root). To use it:
%
%   (a) Launch MATLAB with this repo folder as the startup/current folder, OR
%   (b) run `startup` manually once per session from this folder, OR
%   (c) add this line to your global userpath startup.m so it always loads:
%         run('C:\Users\chany\Documents\gits\Knowledge Enhancement\scFormationRelativeOrbit-CUBoulder\startup.m');
%
% After it runs, the `astro` package is available everywhere:
%       import astro.*
%       [a,e,i,RAAN,omega,f] = RV2COE(r,v,mu);
%   or  [r,v] = astro.COE2RV(a,e,i,RAAN,omega,f,mu);   % no import needed

repoRoot = fileparts(mfilename('fullpath'));
libDir   = fullfile(repoRoot, 'lib');   % parent of +astro must be on the path
addpath(libDir);
fprintf('[startup] Registered shared libraries: %s\n', libDir);
