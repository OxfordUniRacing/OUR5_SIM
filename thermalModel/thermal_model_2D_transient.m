function [T, model] = thermal_model_2D_transient(T_init,t,heat_cell,car_velocity)
    %% 1) Make a transient model
    model = createpde('thermal','transient');
    
    [dl,bt,sf,ig,names] = define_geometry();
    geometryFromEdges(model,dl);
    
    figure('units','normalized','position',[0.15 0.15 0.6 0.6]);
    pdegplot(model,'FaceLabels','on','EdgeLabels','off');
    axis equal
    title('Geometry');


    %% 2) Material properties (you already provide rho & Cp — required for transient)
    conductivity_aluminium = 237;
    density_aluminium = 2700;
    specific_heat_aluminium = 900;
    
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
    [al_faces, ct_faces, cell_faces] = findMaterialFaces(model, ig);
    
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
    

%% 3) Initial conditions (uniform 25 °C — change if you need)
thermalIC(model, T_init);

%% 4) Boundary conditions (same definitions are fine for transient)
[freeEdges, heaterEdges, cellFaces] = getBCindexes(model,ig);

% Convection on underside
htc = flatPlateAirHTC(car_velocity,ig.pack_length);
thermalBC(model,'edge',freeEdges,...
          'ConvectionCoefficient',htc,...
          'AmbientTemperature',25);

% appply heatflux on edges
heat_row_region = heat_cell * ig.n_cell_module / ig.Nrows; % half heatflux as assumed equal comming out either end
Q = heat_row_region * ig.numU * ig.Nrows;
internalHeatSource(model, heat_row_region/(ig.cell_diameter*ig.cell_length), 'face', cellFaces); %heat input is apparently in W/m2 but i am unsure about this and it is applied equally to each region
% --- OR ---
% (B) Time-dependent example: step from 0 to heatflux_cell at t = 60 s
% hf = @(region,state) (state.time >= 60) * heatflux_cell;
% thermalBC(model, 'edge', heaterEdges, 'HeatFlux', hf);

%% 6) Mesh
generateMesh(model,"Hmax",0.5e-3);

%% 7) Time vector and solve
% Choose times you want results at (e.g., 0 to 3600 s every 1 s)

result = solve(model, t);
T = result.Temperature;

