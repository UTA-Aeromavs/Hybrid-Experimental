clc; clear; close all;

% Calculate the motor burn for a real motor
% Tests against table in Humble Page 382 (402 pdf)
a = 2.066e-5;
m = -0.15;
n = 0.75;
rho_s = 1000;
dm_oxdt = 7.95;

N = 100;
Length = 0.5;
di = 0.152;

ri = ones(N, 1) * di / 2;
dx = Length / N;
x = dx:dx:Length;

IB = @(t,r) descrete_internal_ballistics(t, r, dx, dm_oxdt, rho_s, a, m, n);
times = linspace(0, 10, 100);
[times,r] = ode45(IB, times, ri);

% Calculate average radius at each time
r_avg = mean(r, 2);

% Average regression rate
rdot = gradient(r_avg, times);

A = pi*r_avg.^2;
G_ox = dm_oxdt ./ A;

DH_avg = 2*r_avg;

L_bore = Length * ones(size(G_ox));
rho_fg = rho_s * ones(size(G_ox));

scatter(G_ox, rdot)
xlabel('Oxidizer mass flux [kg/m^2/s]')
ylabel('Regression rate [m/s]')
grid on

X = [G_ox DH_avg L_bore rho_fg];

BETA0 = [a m n];

mdl = fitnlm(X, rdot, @mean_internal_ballistics, BETA0)