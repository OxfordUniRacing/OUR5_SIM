function [T, model] = thermal_model_2D_transient(T_init,t,heat_cell,car_velocity)
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

%% 4) Boundary conditions (same definitions are fine for transient)
[freeEdges, heaterEdges] = getBCedges(model,ig);

% Convection on underside
htc = flatPlateAirHTC(car_velocity,ig.pack_length);
thermalBC(model,'edge',freeEdges,...
          'ConvectionCoefficient',htc,...
          'AmbientTemperature',25);

% appply heatflux on edges
heatflux_cell = heat_cell * ig.n_cell_module / ig.Nrows / 2 / ig.heater_height / ig.pack_length; % half heatflux as assumed equal comming out either end
A = ig.heater_height*ig.pack_length*ig.Nrows*2*ig.numU;
Q = A*heatflux_cell;
thermalBC(model, 'edge', heaterEdges, 'HeatFlux', heatflux_cell);
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

