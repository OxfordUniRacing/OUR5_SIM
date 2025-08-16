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
    
    conductivity_cell        = 30;     % example – depends on cell design
    density_cell             = 2500;   % typical 18650 density ~2.5 g/cm3
    specific_heat_cell       = 1000;   % J/kg-K (approx, varies with chemistry)
    
    % find faces
    al_faces   = find(~contains(names,'thermal') & ~contains(names,'cell')); % everything else
    ct_faces   = find(contains(names,'thermal'));
    cell_faces = find(contains(names,'cell'));

    % apply materials
thermalProperties(model, 'Face', al_faces, ...
    'ThermalConductivity', conductivity_aluminium, ...
    'MassDensity', density_aluminium, ...
    'SpecificHeat', specific_heat_aluminium);

thermalProperties(model, 'Face', ct_faces, ...
    'ThermalConductivity', conductivity_TIM, ...
    'MassDensity', density_TIM, ...
    'SpecificHeat', specific_heat_TIM);

thermalProperties(model, 'Face', cell_faces, ...
    'ThermalConductivity', conductivity_cell, ...
    'MassDensity', density_cell, ...
    'SpecificHeat', specific_heat_cell);
    
    %% 3) Initial conditions (uniform 25 °C — change if you need)
    thermalIC(model, T_init);
    
    %% 4) setup boundary conditions (same definitions are fine for transient)
    [freeEdges, heaterEdges, cellFaces] = getBCindexes(model,ig);
    
    [t_unique, ia] = unique(t);
    % Convection on underside
    htc = flatPlateAirHTC(car_velocity,ig.pack_length);
    htc_unique = htc(ia);
    
    % appply heatflux on edges
    heat_row_region = heat_cell * ig.n_cell_module / ig.Nrows; % half heatflux as assumed equal comming out either end
    heat_cellregion_unique = heat_row_region(ia);
    
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
        Qrow_fun  = @(region, state) interp1(t_unique, heat_cellregion_unique, state.time + t_offset, 'linear', 'extrap');
        htc_fun = @(region, state) interp1(t_unique, htc_unique, state.time + t_offset, 'linear', 'extrap');
    
        % Apply BCs
        internalHeatSource(model, Qrow_fun, 'face', cellFaces);
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


