% thermal sim
T_init = params.ambient_temperature;

% average heat and velocity over the lap
for lapN = 1:Num_Laps
    heat_cell(lapN) = mean((I_data(n_steps*(lapN-1) + 1:n_steps + n_steps*(lapN-1))/ params.battery.Np).^2 * params.cellR);
    car_velocity_lapAvg(lapN) = mean(v_data((n_steps*(lapN-1) + 1:n_steps + n_steps*(lapN-1))));
    animation_t(lapN) = t_data(n_steps + n_steps*(lapN-1)); 
end
animation_t = [0,animation_t];
heat_cell_lapAvg = [0,heat_cell_lapAvg];
car_velocity_lapAvg = [0,car_velocity_lapAvg];

% simulation that computes temp profile at end of each lap using an
% averaged heat flux and average car velocity
[T_animations, model] = thermal_model_2D_transient(T_init,animation_t,mean(heat_cell),mean(v_data));

% simulation that uses the profile to update at the model at every time
% point
% [T_animations, model] = thermal_model_2D_transient_profile(T_init,t_data,heat_cell,v_data,animation_t);



%% Plot at a specific time index (e.g., t = 600 s)
k = length(animation_t);   
figure
pdeplot(model,'XYData',T_animations(:,end),'ColorMap','jet');
xlabel("dimensions (m)")
c = colorbar; c.Label.String = 'Temperature (°C)';
[maxT, idx] = max(T_animations(:,end));
title(sprintf('Temperature at t = %.1f s\n Max temperature: %2.2f °C', animation_t(k),maxT));
% add max temp point annotation
fprintf("Max temperature: %2.2f °C\n",maxT);


% Create animated GIF
gifFile = 'thermal_animation.gif';
cbar_interval = 1;
figure
for k = 1:length(animation_t)
    pdeplot(model,'XYData',T_animations(:,k),'ColorMap','jet');
    clim([T_init ceil(max(max(T_animations))/cbar_interval)*cbar_interval])
    c = colorbar; c.Label.String = 'Temperature (°C)';
    [maxT, idx] = max(T_animations(:,k));
    title(sprintf('Temperature at t = %.1f s\n Max temperature: %2.2f °C', animation_t(k),maxT));
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


