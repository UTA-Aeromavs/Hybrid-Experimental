function Engine = combustor_model(ID, L, mdot_ox)
    % Optimizable combustor model
    %     A combustor model that can be optimized for its inputs

    % Get CEA table
    cea = readtable("CEA_Isobaric.csv");
    % Parameters
    p = 2.344e6; % For design purposes we assume this is constant, which is based truly off throat conditons
    OD = 0.038; % This is limited by the physical casing and caps the impulse
    
    % Simplified IB model
    function [value, isterminal, direction] = burn_complete(t, r)
        value = r - 0.038/2;
        isterminal = 1;
        direction = 0;
    end

    a = 0.188;
    n = 0.347;
    G_ox = @(r) mdot_ox / pi / r^2;
    rdot = @(t, r) simplified_internal_ballistics([a, n], G_ox(r))/1000;
    options = odeset('Events', @burn_complete);
    [t, r] = ode45(rdot, [0 10], ID/2, options);

    % Calculate performance vectors
    rho = 920;
    S = arrayfun(@(r) pi*r*2*L, r);
    A = arrayfun(@(r) pi*r^2, r);
    mdot_f = arrayfun(@(r, S) rdot(0, r)*S*rho, r, S);
    of = mdot_ox ./ mdot_f;
    mdot = mdot_ox+mdot_f;
    F = arrayfun(@(of, mdot) mdot * interp1(cea.OF, cea.ve, of), of, mdot);

    % Calculate scalar performance characteristics
    It = trapz(t, F);
    of_avg = trapz(t, of)/t(end);
    m_f = rho*pi*(OD^2 - ID^2)/4*L;
    m_ox = mdot_ox * t(end);
    m = m_f + m_ox;
    Is = It / m;

    % Output
    Engine.ID = ID;
    Engine.OD = OD;
    Engine.L = L;
    Engine.p = p;
    Engine.rV = r;
    Engine.GoxV = mdot_ox ./ A;
    Engine.GV = mdot ./ A;
    Engine.SV = S;
    Engine.AV = A;
    Engine.tV = t;
    Engine.mdot_fV = mdot_f;
    Engine.ofV = of;
    Engine.It = It;
    Engine.of = of_avg;
    Engine.m_f = m_f;
    Engine.m_ox = m_ox;
    Engine.m = m;
    Engine.Is = Is;
end
