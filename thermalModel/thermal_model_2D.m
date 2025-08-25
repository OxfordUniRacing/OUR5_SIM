function [T, model] = thermal_model_2D(T_init,t,heat_cell,car_velocity,sim_type)
    %% 1) Make a transient model
    if strcmp(sim_type,'transient')
        model = createpde('thermal','transient');
    elseif strcmp(sim_type,'steadystate')
        model = createpde('thermal','steadystate');
    else
        error("invalid sim_type, transient or steadystate")
    end

    [dl,bt,sf,ig,names] = define_geometry_2(true); % true false for fins
    geometryFromEdges(model,dl);
    
    figure('units','normalized','position',[0.15 0.15 0.6 0.6]);
    pdegplot(model,'FaceLabels','on','EdgeLabels','on');
    axis equal
    title('Geometry');


    %% 2) Material properties (you already provide rho & Cp — required for transient)
    conductivity_aluminium = 237;
    density_aluminium = 2700;
    specific_heat_aluminium = 900;

    conductivity_copper = 401;     % W/(m·K)
    density_copper = 8960;         % kg/m^3
    specific_heat_copper = 385;    % J/(kg·K)
    copper_region_scaling_factor = (ig.NbusbarModule/2 * ig.busbar_width) / ig.pack_length;
    
    conductivity_TIM = 3;
    density_TIM = 3000;
    specific_heat_TIM = 100;
    
    cell_volume = (pi * ig.cell_diameter^2/4) * ig.cell_length;
    cell_region_volume = ig.cell_diameter * ig.cell_length * ig.pack_length;
    conductivity_cell        = 10;     % example – depends on cell design
    density_cell             = 0.07 / cell_volume;   % cells weight 70g accoroding to molicel datasheet
    specific_heat_cell       = 830;   % J/kg-K (approx, varies with chemistry, 1400 for NMC, 830 for NCA)
    
    region_scaling_factor    = (cell_volume * ig.n_cell_module / ig.Nrows) / cell_region_volume; % the region assigned to cells is cuboidic but the cells dont fill the full volume obviosuly
    
    conductivity_cellregion  =  conductivity_cell * region_scaling_factor;
    density_cellregion       =  density_cell * region_scaling_factor;
    
    % find faces
    [al_faces, ct_faces, cell_faces, cu_faces] = findMaterialFaces(model, ig);
    
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
    'ThermalConductivity', conductivity_cellregion, ...
    'MassDensity', density_cellregion, ...
    'SpecificHeat', specific_heat_cell);

    if cu_faces
        thermalProperties(model, 'Face', cu_faces, ...
            'ThermalConductivity', conductivity_copper * copper_region_scaling_factor, ...
            'MassDensity', density_copper, ...
            'SpecificHeat', specific_heat_copper * copper_region_scaling_factor);
    end
            

%% 3) Initial conditions (uniform 25 °C — change if you need)
thermalIC(model, T_init);

%% 4) Boundary conditions (same definitions are fine for transient)
[freeEdges, cellFaces] = getBCindexes(model,ig);
htc = flatPlateAirHTC(car_velocity,ig.pack_length);
if strcmp(sim_type, 'transient')
    if isscalar(heat_cell) && isscalar(car_velocity) 
        % Convection on underside
        thermalBC(model,'edge',freeEdges,...
              'ConvectionCoefficient',htc,...
              'AmbientTemperature',25);
    
        % appply heatflux on edges
        heat_row_region = heat_cell * ig.n_cell_module / ig.Nrows; % half heatflux as assumed equal comming out either end
        Q = heat_row_region * ig.numU * ig.Nrows;
        internalHeatSource(model, heat_row_region/(ig.cell_diameter*ig.cell_length)/ig.pack_length, 'face', cellFaces); %heat input is apparently in W/m2 but i am unsure about this and it is applied equally to each region
    
    else
        % Create boundary condition functions using interpolation
        heat_row_region = heat_cell * ig.n_cell_module / ig.Nrows; % half heatflux as assumed equal comming out either end
        Q_fun = @(region, state) ...
        repmat( interp1(t, heat_row_region/(ig.cell_diameter*ig.cell_length)/ig.pack_length, state.time, 'linear', 'extrap'), ...
                1, numel(region.x) );
        htc_fun = @(region, tsim) interp1(t, htc, tsim.time, 'linear', 'extrap');
        
        % Apply BCs
        internalHeatSource(model, Q_fun, 'face', cellFaces);
        thermalBC(model, 'edge', freeEdges, ...
            'ConvectionCoefficient', htc_fun, ...
            'AmbientTemperature', T_init);
    end
else
    % Convection on underside
    thermalBC(model,'edge',freeEdges,...
          'ConvectionCoefficient',htc(end),...
          'AmbientTemperature',T_init);

    % appply heatflux on edges
    heat_row_region = heat_cell(end) * ig.n_cell_module / ig.Nrows; 
    Qtot = heat_cell(end) * ig.numU * ig.n_cell_module;
    internalHeatSource(model, heat_row_region/(ig.cell_diameter*ig.cell_length)/ig.pack_length, 'face', cellFaces); %heat input is apparently in W/m2 but i am unsure about this and it is applied equally to each region
end

%% 6) Mesh
generateMesh(model,"Hmax",0.5e-3);

%% 7) Time vector and solve
% Choose times you want results at (e.g., 0 to 3600 s every 1 s)
if strcmp(sim_type, 'transient')
    result = solve(model, t);
else
    result = solve(model);
end
% evaluateHeatRate(result,"edge",freeEdges)*ig.pack_length can be used to check the output 
% heat matches the generated heat Qtot when last checked it did

T = result.Temperature;

