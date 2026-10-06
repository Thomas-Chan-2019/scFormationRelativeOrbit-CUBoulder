clear all; close all;
import astro.*   % RV2COE / COE2RV from lib/+astro (run `startup` if not found)
% test_COE_RV_roundtrip
% Verifies COE2RV and RV2COE are mutual inverses.
%
% Two independent checks:
%   (A) COE -> RV -> COE : recovered elements match the originals.
%   (B) RV  -> COE -> RV : recovered state matches the original vectors.
%
% A clean round-trip in BOTH directions confirms the leading -mu/h sign in
% COE2RV's velocity is internally consistent with RV2COE. A velocity-only
% failure (position OK, velocity flipped/wrong) would point at that sign.

mu = 3.986e14;              % [m^3/s^2] Earth
tol = 1e-6;                % relative tolerance

%% (A) COE -> RV -> COE
% Reference elements (Quiz 10 Q6 style, propagated to a true anomaly).
a0     = 8000e3;           % [m]
e0     = 0.1;
i0     = 30*pi/180;        % [rad]
RAAN0  = 145*pi/180;       % [rad]
omega0 = 120*pi/180;       % [rad]
f0     = 47.35*pi/180;     % [rad] arbitrary true anomaly

[r_vec, v_vec] = COE2RV(a0,e0,i0,RAAN0,omega0,f0,mu);
[a,e,i,RAAN,omega,f] = RV2COE(r_vec,v_vec,mu);

% Wrap angle differences into [-pi,pi] before comparing.
wrap = @(d) mod(d + pi, 2*pi) - pi;

errA = [ (a-a0)/a0 ; ...
         (e-e0) ; ...
         wrap(i-i0) ; ...
         wrap(RAAN-RAAN0) ; ...
         wrap(omega-omega0) ; ...
         wrap(f-f0) ];

fprintf('=== (A) COE -> RV -> COE ===\n');
fprintf('  d a/a0     = % .3e\n', errA(1));
fprintf('  d e        = % .3e\n', errA(2));
fprintf('  d i    [rad]= % .3e\n', errA(3));
fprintf('  d RAAN [rad]= % .3e\n', errA(4));
fprintf('  d omega[rad]= % .3e\n', errA(5));
fprintf('  d f    [rad]= % .3e\n', errA(6));
assert(max(abs(errA)) <= tol, 'Round-trip (A) FAILED: max err %.3e', max(abs(errA)));
fprintf('  PASS (max err %.3e)\n\n', max(abs(errA)));

%% (B) RV -> COE -> RV
% Reference state (Quiz 10 Q7).
r0 = [-820.865, -1905.95, -7445.9]' * 1e3;   % [m]
v0 = [-6.75764, -1.85916, 0.930651]' * 1e3;  % [m/s]

[a,e,i,RAAN,omega,f] = RV2COE(r0,v0,mu);
[r_vec, v_vec] = COE2RV(a,e,i,RAAN,omega,f,mu);

err_r = norm(r_vec - r0)/norm(r0);
err_v = norm(v_vec - v0)/norm(v0);

fprintf('=== (B) RV -> COE -> RV ===\n');
fprintf('  ||dr||/||r0|| = % .3e\n', err_r);
fprintf('  ||dv||/||v0|| = % .3e\n', err_v);
assert(err_r <= tol, 'Round-trip (B) position FAILED: %.3e', err_r);
assert(err_v <= tol, 'Round-trip (B) velocity FAILED: %.3e  <-- check -mu/h sign', err_v);
fprintf('  PASS (dr %.3e, dv %.3e)\n\n', err_r, err_v);

fprintf('ALL ROUND-TRIP CHECKS PASSED.\n');
