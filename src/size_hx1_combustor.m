clc; clear; close all;

% Optimization Variables
ID = optimvar('ID', ...
    LowerBound=0.001, UpperBound=0.038);       % [m]

L = optimvar('L', ...
    LowerBound=0.001, UpperBound=0.140);       % [m]

mdot_ox = optimvar('mdot_ox', ...
    LowerBound=0.001);                         % [kg/s]

% Model Performance
It = fcn2optimexpr( ...
    @(ID,L,mdot_ox) getfield( ...
        combustor_model(ID,L,mdot_ox), 'It'), ...
    ID, L, mdot_ox, Analysis='off');

of_variation = fcn2optimexpr( ...
    @of_std, ID, L, mdot_ox, Analysis='off');

% Optimization Problem
prob = optimproblem(ObjectiveSense='minimize');

prob.Objective = of_variation;

% Maximum allowable impulse [N*s]
prob.Constraints.impulse = It == 2000;

% Maximum oxidizer mass flux [kg/(m^2*s)]
prob.Constraints.G_ox = 4*mdot_ox <= 600*pi*ID^2;

% Initial Guess
x0.ID = 0.022;
x0.L = 0.138;
x0.mdot_ox = 0.092;

% Solver Options
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'PlotFcn', @optimplotfvalconstr);

% Solve
tic
[sol, fval, exitflag, output] = solve(prob, x0, Options=options);
toc

% Evaluate Optimized Design
Engine = combustor_model(sol.ID, sol.L, sol.mdot_ox);

disp(sol)
disp(Engine)

% Objective Function
function sigma = of_std(ID, L, mdot_ox)

    Engine = combustor_model(ID, L, mdot_ox);

    t = Engine.tV;
    of = Engine.ofV;

    of_avg = trapz(t, of)/t(end);

    sigma = sqrt(trapz(t, (of - of_avg).^2)/t(end));

end
