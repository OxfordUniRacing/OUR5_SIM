function h = flatPlateAirHTC(v, L)
    % flatPlateAirHTC  Calculate average heat transfer coefficient for air over flat plate
    % v [m/s]  - flow velocity
    % L [m]    - plate length
    % h [W/m^2/K] - average heat transfer coefficient

    % --- Air properties at ~25°C ---
    rho = 1.184;              % kg/m^3
    nu = 15.06e-6;             % m^2/s
    k = 0.0257;                % W/m/K
    Pr = 0.71;                 % Prandtl number

    % --- Reynolds number ---
    ReL = v * L / nu;

    % --- Nusselt number ---
    if ReL <= 5e5
        % Laminar
        NuL = 0.664 .* ReL.^0.5 .* Pr.^(1/3);
    else
        % Turbulent
        NuL = 0.037 .* ReL.^(4/5) .* Pr.^(1/3);
    end

    % --- Heat transfer coefficient ---
    h = NuL * k / L;
end