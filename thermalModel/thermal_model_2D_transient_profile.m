function [T_animations, model] = thermal_model_2D_transient_profile(T_init,t,heat_cell,car_velocity,t_animation)
    %% 1) Make the model transient (instead of 'steadystate')
    model = createpde('thermal','transient');
    
    [dl,bt,sf,ig,names] = define_geometry();
    geometryFromEdges(model,dl);
    
    figure('units','normalized','position',[0.15 0.15 0.6 0.6]);
    pdegplot(model,'FaceLabels','on','EdgeLabels','on');
    axis equal
    title('Geometry');
    
    %% 2) Material properties (you already provide rho & Cp — required for transient)
    condutivity_aluminium = 237;
    density_aluminium = 2700;
    specific_heat_alumiunium = 900;
    
    condutivity_TIM = 10;
    density_TIM = 500;
    specific_heat_IIM = 100;
    
    %find faces
    al_faces = find(~contains(names,'thermal'));
    ct_faces = find(contains(names,'thermal'));
    
    %apply materials
    thermalProperties(model, 'Face', al_faces, ...
        'ThermalConductivity', condutivity_aluminium, ...
        'MassDensity', density_aluminium, ...
        'SpecificHeat', specific_heat_alumiunium);
    
    thermalProperties(model, 'Face', ct_faces, ...
        'ThermalConductivity', condutivity_TIM, ...
        'MassDensity', density_TIM, ...
        'SpecificHeat', specific_heat_IIM);
    
    
    %% 3) Initial conditions (uniform 25 °C — change if you need)
    thermalIC(model, T_init);
    
    %% 4) setup boundary conditions (same definitions are fine for transient)
    [freeEdges, heaterEdges] = getBCedges(model,ig);
    
    [t_unique, ia] = unique(t);
    % Convection on underside
    htc = flatPlateAirHTC(car_velocity,ig.pack_length);
    htc_unique = htc(ia);
    
    % appply heatflux on edges
    heatflux_cell = heat_cell * ig.n_cell_module / ig.Nrows / 2 / ig.heater_height / ig.pack_length; % half heatflux as assumed equal comming out either end
    heatflux_cell_unique = heatflux_cell(ia);
    
    %% 6) Mesh
    generateMesh(model,"Hmax",0.5e-3);
    
    %% 7) Time vector and solve     
    % --- animation loop ---
    for i_save = 1:length(t_animation)-1
        % Save results here
        % Local time vector for each lap
        t_chunk = t_unique(find(t_unique>=t_animation(i_save)):find(t_unique>=t_animation(i_save+1)));
        % Global time offset for this lap
        t_offset = t_chunk(1);
    
        % Create boundary condition functions using interpolation
        hf_fun  = @(region, state) interp1(t_unique, heatflux_cell_unique, state.time + t_offset, 'linear', 'extrap');
        htc_fun = @(region, state) interp1(t_unique, htc_unique, state.time + t_offset, 'linear', 'extrap');
    
        % Apply BCs
        thermalBC(model, 'edge', heaterEdges, 'HeatFlux', hf_fun);
        thermalBC(model, 'edge', freeEdges, ...
            'ConvectionCoefficient', htc_fun, ...
            'AmbientTemperature', T_init);
    
        % Solve this lap
        result = solve(model, t_chunk);
    
        % Save only final state
        T_animations(:, i_save) = result.Temperature(:, end);
    
        % Set final state as initial for next lap
        nodes = result.Mesh.Nodes;                 % 2 x N in 2-D
        F = scatteredInterpolant(nodes(1,:).', nodes(2,:).', ...
                                  result.Temperature(:,end), ...
                                  'linear','nearest');
    
        % Set spatially varying IC for the next chunk
        thermalIC(model, @(location) F(location.x, location.y));

    end


end


