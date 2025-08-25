% thermal sim
T_init = params.ambient_temperature;

% average heat and velocity over the lap
clear heat_cell_lapAvg car_velocity_lapAvg lap_t
for lapN = 1:Num_Laps
    heat_cell_lapAvg(lapN) = mean((I_data(n_steps*(lapN-1) + 1:n_steps + n_steps*(lapN-1))/ params.battery.Np).^2 * params.cellR);
    car_velocity_lapAvg(lapN) = mean(v_data((n_steps*(lapN-1) + 1:n_steps + n_steps*(lapN-1))));
    lap_t(lapN) = t_data(n_steps + n_steps*(lapN-1)); 
end
lap_t = [0,lap_t];
heat_cell_lapAvg = [0,heat_cell_lapAvg];
car_velocity_lapAvg = [0,car_velocity_lapAvg];

% simulation that computes temp profile at end of each lap using an
% averaged heat flux and average car velocity
[lap_T, model] = thermal_model_2D(T_init,lap_t,heat_cell_lapAvg,car_velocity_lapAvg,'transient');

% simulation that uses the profile to update at the model at every time
% point
% [T_animations, model] = thermal_model_2D_transient_profile(T_init,t_data,heat_cell,v_data,animation_t);



%% Plot at a specific time index (e.g., t = 600 s)
k = length(lap_t);   
figure
pdeplot(model,'XYData',lap_T(:,end),'ColorMap','jet');
xlabel("dimensions (m)")
c = colorbar; c.Label.String = 'Temperature (°C)';
[maxT, idx] = max(lap_T(:,end));
title(sprintf('Temperature at t = %.1f s\n Max temperature: %2.2f °C', lap_t(k),maxT));
% add max temp point annotation
fprintf("Max temperature: %2.2f °C\n",maxT);


% Create animated GIF
gifFile = 'thermal_animation.gif';
cbar_interval = 1;
figure
for k = 1:length(lap_t)
    pdeplot(model,'XYData',lap_T(:,k),'ColorMap','jet');
    clim([T_init ceil(max(max(lap_T))/cbar_interval)*cbar_interval])
    c = colorbar; c.Label.String = 'Temperature (°C)';
    [maxT, idx] = max(lap_T(:,k));
    title(sprintf('Temperature at t = %.1f s\n Max temperature: %2.2f °C', lap_t(k),maxT));
    axis equal
    drawnow

    % Capture the frame as an image
    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);

    if k == 1
        imwrite(A,map,gifFile,'gif','LoopCount',Inf,'DelayTime',0.5);
    else
        imwrite(A,map,gifFile,'gif','WriteMode','append','DelayTime',0.5);
    end
end


%% hot spot plot
[hotspotEndTemp, hotspotIdx] = max(lap_T(:,end));   % 1xN vector (max over nodes)
hotspotTemp = lap_T(hotspotIdx,:);
figure
plot(lap_t, hotspotTemp, 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Hot spot temperature (°C)')
grid on
title(sprintf('Maximum (Hot Spot) Temperature vs Time\n Final Temperature = %2.2f °C',hotspotEndTemp))

%% thermal resistance
[steadystate_T, model] = thermal_model_2D(T_init,lap_t,heat_cell_lapAvg,car_velocity_lapAvg,'steadystate');
figure
pdeplot(model,'XYData',steadystate_T,'ColorMap','jet');
xlabel("dimensions (m)")
c = colorbar; c.Label.String = 'Temperature (°C)';
[maxT, idx] = max(steadystate_T);
title(sprintf('SteadyState, Max temperature: %2.2f °C',maxT));
hotspotEndTemp = max(steadystate_T(:));   % 1xN vector (max over nodes)
maxRth = (hotspotEndTemp-T_init)/heat_cell_lapAvg(end)