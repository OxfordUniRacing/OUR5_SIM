function [result, model] = thermal_model_2D_transient_profile(T_init,t,heat_cell,car_velocity,t_animation)

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
[t_unique, ia] = unique(t);
htc_unique = htc(ia);
htc_fun = @(region,state) interp1(t_unique, htc_unique, state.time, 'linear', 'extrap');
thermalBC(model, 'edge', freeEdges, ...
          'ConvectionCoefficient', htc_fun, ...
          'AmbientTemperature', T_init);

% appply heatflux on edges
heatflux_cell = heat_cell * ig.n_cell_module / ig.Nrows / 2 / ig.heater_height / ig.pack_length; % half heatflux as assumed equal comming out either end
[t_unique, ia] = unique(t);
heatflux_cell_unique = htc(ia);
% (B) Time-dependent example: step from 0 to heatflux_cell at t = 60 s
hf = @(region,state) interp1(t_unique, heatflux_cell_unique, state.time, 'linear', 'extrap');
thermalBC(model, 'edge', heaterEdges, 'HeatFlux', hf);

%% 6) Mesh
generateMesh(model,"Hmax",0.5e-3);

%% 7) Time vector and solve
% Choose times you want results at (e.g., 0 to 3600 s every 1 s)
result = solve(model, unique(t));

%% 8) Plot at a specific time index (e.g., t = 600 s)
% k = length(t);   % or k = 601 if you know the index
% figure
% pdeplot(model,'XYData',result.Temperature(:,k),'ColorMap','jet');
% xlabel("dimensions (m)")
% c = colorbar; c.Label.String = 'Temperature (°C)';
% title(sprintf('Temperature at t = %.1f s', t(k)));
% 
% %% 9) Create animated GIF
% gifFile = 'thermal_animation.gif';
% 
% for k = 1:length(t_animation)
%     pdeplot(model,'XYData',result.Temperature(:,k),'ColorMap','jet');
%     caxis([T_init ceil(max(max(result.Temperature))/10)*10])
%     c = colorbar; c.Label.String = 'Temperature (°C)';
%     title(sprintf('Temperature at t = %.2f s',t(k)));
%     drawnow
% 
%     % Capture the frame as an image
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [A,map] = rgb2ind(im,256);
% 
%     if k == 1
%         imwrite(A,map,gifFile,'gif','LoopCount',Inf,'DelayTime',0.5);
%     else
%         imwrite(A,map,gifFile,'gif','WriteMode','append','DelayTime',0.5);
%     end
% end
end