function drdt = mean_internal_ballistics(beta, x)
    % Mean Internal Ballistics Model
    %   Determines the regression rate using Humble Average internal Ballistics Model
    %   Can be used in fitnlm to find the internal ballistics Coefficient estimate given average IB
    %   values
    %
    % Inputs
    %   beta    array of regression rate coefficients
    %   (1) a       IB Coefficient
    %   (2) m       IB Coefficient
    %   (3) n       IB Coefficient
    %   x       Input Matrix
    %   (1) G_ox    Oxidizer Mass Flux
    %   (2) DH_avg  Hydraulic Diameter
    %   (3) L_bore  Length of Combustion Chamber
    %   (4) rho_fg  Density of fuel Grain
    %
    % Outputs
    %   drdt    Fuel Regression Rates

    % Regression Rate Coefficients
    a = beta(1);
    m = beta(2);
    n = beta(3);
    
    G_ox = x(:, 1); % Oxidizer Mass flux rate
    DH_avg = x(:, 2); % Average Hydrolic Diameter
    L_bore = x(:, 3); % Total Chamber Length
    rho_fg = x(:, 4); % Fuel Grain Density

    % Humble Equation 7.34 
    % Page 384 (404 on PDF)
    drdt = a/(1+m) .* G_ox.^n .* L_bore.^m .*(1+2*n*a/(1+m).* rho_fg.*L_bore.^(1+m)./ DH_avg ./ G_ox.^(1-n));
end