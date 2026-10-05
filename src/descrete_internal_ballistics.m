function [drdt, dmdt, of] = descrete_internal_ballistics(t, r, dx, dm_oxdt, rho_s, a, m, n, r_ch)
    % Descrete Internal Ballistics Model
    %   Calculates the internal ballistics of a hybrid rocket engine
    %
    % Inputs
    %   t       Time (Unused currently)             [s]
    %   r       Evenly spaced array of cell radii   [m]
    %   dx      Width of cells                      [m]
    %   dm_oxdt Oxidizer mass Flow rate             [kg/m^2/s]
    %   rho_s   Density of solid propellant grain   [kg/m^3]
    %   a       IB Coefficient
    %   m       IB Coefficient
    %   n       IB Coefficient
    %   r_ch    Physical radius of combustion chamber
    %
    % Outputs
    %   drdt    Regression rate
    %   dmdt    Combined (oxidizer and fuel) Mass flow rate
    %   of      Oxidizer to fuel mass ratio

    arguments
        t % currently Unused but planned to be used later
        r (:,1) double {mustBePositive}
        dx (:,1) double {mustBePositive}
        dm_oxdt (1,1) double {mustBePositive}
        rho_s (1,1) double {mustBePositive}
        a (1,1) double {mustBePositive}
        m (1,1) double
        n (1,1) double
        r_ch (1,1) double {mustBeGreaterThan(r_ch, r)} = inf
    end

    % Sanity checks
    len = length(r);
    if isscalar(dx)
        dx = repmat(dx, len, 1);
    elseif length(dx) ~= len
        error("Supplied cell dimension is different from radius dimension");
    end
    if any(dx <= 0) || any(r <= 0)
        error("Radius and cell dimensions must be positive");
    end
    drdt = zeros(len, 1);

    % Oxidizer mass flow entering grain
    dmdt = dm_oxdt;
    x=0;

    for j = 1:len
        x = dx(j) + x;

        % Local Parameters
        A = pi*r(j)^2;
        P = 2*pi*r(j);
        G = dmdt/A;

        % Local Internal Ballistics
        drdt(j) = a * G^n * x^m;
        if r(j) - drdt(j) > r_ch
            drdt(j) = r(j) - r_ch;
        end

        % Fuel generated in this axial cell
        dmdot_f = rho_s * P * drdt(j) * dx(j);

        % Carry total mass flow to next cell
        dmdt = dmdt + dmdot_f;
    end
    of = dm_oxdt / (dmdt - dm_oxdt);
end

function out = set_size(in, other)
    if isscalar(in)
        out = repmat(in, size(other));
    elseif size(in) == size(other)
        out = in;
    else
        error('Size different between output and input is undefined');
    end
end