car_v = 15;
heat_cell = 3

%%%%% convert to geometry and import to finite element model
model = createpde('thermal','steadystate');
[dl,bt,sf,ig,names] = define_geometry();
geometryFromEdges(model,dl);

figure('units','normalized','position',[0.15 0.15 0.6 0.6]);
pdegplot(model,'FaceLabels','on','EdgeLabels','on');
axis equal
title('Geometry');


%%%%% material properties
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

%%%%% APPPLY BOUNDARY CONDITIONS
[freeEdges, heaterEdges] = getBCedges(model,ig);

% Convection on underside
htc = flatPlateAirHTC(car_v,ig.pack_length);
thermalBC(model,'edge',freeEdges,...
          'ConvectionCoefficient',htc,...
          'AmbientTemperature',25);

% appply heatflux on edges
heatflux_cell = heat_cell * ig.n_cell_module / ig.Nrows / 2 / ig.heater_height / ig.pack_length; % half heatflux as assumed equal comming out either end
A = ig.heater_height*ig.pack_length*ig.Nrows*2*ig.numU;
Q = A*heatflux_cell
thermalBC(model, 'edge', heaterEdges, 'HeatFlux', heatflux_cell);

%%%% Solve and plots
generateMesh(model,"Hmax",0.5e-3)
result = solve(model);

figure
pdeplot(model,'XYData',result.Temperature,'ColorMap','jet');
xlabel("dimensions (m)")
% axis equal 
c = colorbar;
c.Label.String = 'Temperature (°C)';
