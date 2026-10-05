function drdt = simplified_internal_ballistics(beta, G_ox)
    % Simplified Internal Ballisitics Model
    %   Calculates the simple internal ballistics model found in most papers aG_ox^n
    %
    % Inputs
    %   beta    array of regression rate coefficients
    %   (1) a       IB Coefficient
    %   (2) n       IB Coefficient
    %
    % Outputs
    %   G_ox    Oxidizer mass flux

    a = beta(1);
    n = beta(2);
    drdt = a .* G_ox.^n;
end