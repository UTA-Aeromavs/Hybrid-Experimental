clc; clear; close all;
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
times = [0.1 20 40 60];
[times,r] = ode45(IB, times, ri);

internal_diameter_vs_length = figure;
hold on;
xlabel('Distance Down Port x [m]');
ylabel('Port Diameter [m]');
legend('show');

regression_rate_vs_length = figure;
hold on;
xlabel('Distance Down Port x [m]');
ylabel('Regression Rate [m/s]');
legend('show');

for idx = 1:size(r, 1)
    row = r(idx, :);
    t = times(idx);
    [drdt, dmdt(idx), of(idx)] = descrete_internal_ballistics(0, row, dx, dm_oxdt, rho_s, a, m, n);
    figure(internal_diameter_vs_length);
    hold on;
    plot(x, row.*2, 'DisplayName',['t = ', num2str(t), 's']);
    figure(regression_rate_vs_length);
    hold on;
    plot(x, drdt, 'DisplayName',['t = ', num2str(t), 's']);
end

mass_flow_rate_vs_time = figure;
hold on;
xlabel('Time');
ylabel('mass flow rate');
plot(times, dmdt);

of_vs_time = figure;
hold on;
xlabel('Time');
ylabel('Oxidizer to fuel ratio');
plot(times, of);
