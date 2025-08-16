%load_tranfer_model.m
%Basic point mass sim for determining battery spec
%Lewis Blake

clearvars
close all

load('curve.mat','curv', 'dels','track_length');

addpath("./thermalModel")

%CAR PARAMETERS
params.M = 320; %kg
params.M_dist = 0.5; %distribution mass over front wheels
params.gratio = 4.5; %gear reduction ratio
params.lat_mu = 1.5; %lateral tyre coeff. friction
params.long_mu = 1.3; %longitudinal tyre coeff. friction
params.tyre_dia = 16; %tyre diameter, inches
params.COG_h = 0.3; %m
params.wheelbase = 1.525; %m

%aero
params.Cd= 1; % drag co-efficient
params.Cl = 0.1;% lift co-efficient (NOT IMPLEMENTED)
params.frontal_area = 1.2; %m^2
params.air_density = 1.225; 

% BATTERY PARAMETERS
params.max_charge_Crate = 1; % max charging C-rate
params.max_discharge_Crate = 10; % max charging C-rate 
params.battery.Np = 5; % number of cells in parallel in a cell group
params.cell_Ah = 4.5; % battery pack Amp-hours
params.pack_Ah = params.cell_Ah * params.battery.Np; % battery pack Amp-hours
params.cellR = 15e-3; % cell resistance 
params.cellV = 4.2;
params.battery.Ns = 90;
params.battery.cell_specific_heat = 830; % 830 typical for a NCA cell, 1040 would be typical for NMC
params.battery.cell_mass = 70e-3; % from molicell datasheet, 70g
params.battery.cellRth = 10;
%EFFICIENCIES
%Motor efficiency is in seperate "motor_efficiecy.m"
params.efficiency.mechanical = 0.92;

%CONTROL
params.control.driver_skill = 0.8; %Driver skill factor (~0.5 to 1), acts as derate
params.control.driver_smoothness_alpha = 0.95; % smoothing factor, 0 = slow change, 1 = instant change
params.control.max_power = 80e3;
params.control.regen = true; % toggle regen on or off

%ENVIRONMENT
g = 9.81;
params.ambient_temperature = 25;

%TRACK PARAMETERS
endurance_length = 22e3;
num_laps_endurance = 27;
lap_length = endurance_length/num_laps_endurance; %m
corner_min_rad = 5; %m

Num_Laps = 27; % number of simulated laps (27 laps in endurance)


%TRACK SCALING
curv_scale_val = (1/corner_min_rad) / abs(max(abs(curv))); %min hairpin radius of Xm
curv_scale = curv * curv_scale_val;
length_scale = lap_length / track_length;
dels_scale = dels * length_scale;

%INITIALISATION
state.v = 5;
state.Fz = params.M*g;
state.Fz_drive = params.M*g*(1-params.M_dist);
state.F = 0;
state.F_long_load_transfer = 0;
state.a_long = 0;
state.SoC = 1; % Pack state of charge as a fraction
state.battery_voltage = pack_voltage(params,state);
state.brake_flag = 0;
state.brake_index = 0;
state.cell_temperature = params.ambient_temperature;

%storage(1) = state;

%DERIVED
Fy_max = state.Fz*params.lat_mu;

%SIM SETUP
sim.runover = 50;

%MAX VELOCITY (FOR COMPARISON)
x = 1:length(curv_scale);

for k = x
    velocity_max(k) = sqrt(Fy_max/(abs(curv_scale(k)) * params.M));
end

TF = islocalmin(velocity_max);
%plot(x,velocity_max,x(TF),velocity_max(TF),'r*')

%NECESSARY FUNCTIONS
function F_lateral = cornering(state,curv,params)
    F_lateral = params.M * state.v^2 * abs(curv);
end

function F_brake = brake(F_lateral,state,params)
    F_brake = state.Fz * params.long_mu * sqrt(1-(F_lateral/(state.Fz*params.lat_mu)));
    
    if abs(imag(F_brake)) > 0
        F_brake = 0;
    end
    
    F_brake = params.control.driver_skill * F_brake; %  get the max braking force
    F_brake = (1 - params.control.driver_smoothness_alpha) * -state.F + params.control.driver_smoothness_alpha * F_brake;% Smooth braking force (F_brake is positive state.F is negative for braking)

end

function RPM_motor = motor_rpm(state,params)
    RPM_motor = 60 * state.v * params.gratio / (params.tyre_dia * 25.4 * 10^-3 * pi);
end

function [F_drive, T_motor, efficiency, grip_limited] = drive(F_lateral,state,params)
    F_drive_grip = state.Fz_drive * params.long_mu * sqrt(1-(F_lateral/(state.Fz*params.lat_mu)));
    if abs(imag(F_drive_grip)) > 0
        F_drive_grip = 0;
    end
    RPM_motor = motor_rpm(state,params);

    T_motor_max = max_torque(RPM_motor, params,state);
    T_wheel_max = T_motor_max * params.efficiency.mechanical  * params.gratio;
    F_wheel_max = T_wheel_max / (params.tyre_dia * 25.4 * 10^-3 / 2);

    if F_drive_grip < F_wheel_max
        F_drive = F_drive_grip;
        F_drive = F_drive * params.control.driver_skill; %apply drive skill factor if grip limited (not required for getting to max torque limit)
        grip_limited = 1;
    else
        F_drive = F_wheel_max;
        grip_limited = 0;
    end
    % smooth the driving force
    F_drive = (1 - params.control.driver_smoothness_alpha) * state.F + params.control.driver_smoothness_alpha * F_drive;% Smooth driving force
    % compute torque
    T_motor = F_drive * (params.tyre_dia * 25.4 * 10^-3 / 2) / (params.gratio * params.efficiency.mechanical);
    efficiency = motor_efficiency(RPM_motor,T_motor);
end

function [F_long_load_transfer, a_long] = long_load_transfer(params,storage)
    a_long = storage(end).F / params.M;
    F_long_load_transfer = params.M * a_long * params.COG_h / params.wheelbase;
end

%SIM
for lapN = 1:Num_Laps
    for i = 1:length(curv_scale)

        % make sure that the driver doesnt come off the brake before the
        % corner apex which would cause oscillation
        if state.brake_flag == 1 && i < state.brake_index 
            state.brake_flag = 1;
        else
            state.brake_flag = 0;
        end

        state.t = dels_scale(i) / state.v;
        state.Fz_drive = params.M*g*(1-params.M_dist) + state.F_long_load_transfer;
        F_lateral = cornering(state,curv_scale(i),params);
        [F_drive, T_motor, efficiency, grip_limited] = drive(F_lateral,state,params);
        F_aero = aero_force(state,params);

        T_motor = T_motor * params.control.driver_skill;

        F_vehicle = F_drive - F_aero;
        v_trial = state.v + state.t * F_vehicle / params.M;

        state.t = dels_scale(i) / state.v; % time taken to travel through the track segment
        state.Fz_drive = params.M * g * (1-params.M_dist) + state.F_long_load_transfer; % vertical load on the driven wheels
        F_lateral = cornering(state, curv_scale(i), params); % lateral load on the tyres
        [F_drive, T_motor, Eff_motor, grip_limited] = drive(F_lateral,state,params); % drive function returns key variables for ideal drive 
   
        F_vehicle = F_drive - F_aero;

        if state.brake_flag == 0 % run look ahead braking check if not already breaking
            temp = state;
            tempPrev = temp;
            temp.v = v_trial;
            % iterates through a forward looking lap of the track from the current segment plus an
            % additional runover distance to check if the driver will lose control
            % in a corner at the current speed. checks if the driver needs to start
            % braking
        
            for k = i+1:length(curv_scale)+sim.runover
                % if index exceeds track length loop back to the start of the track
                if k > length(curv_scale)
                    o = k - length(curv_scale);
                else
                    o = k;
                end
                
                temp.RPM_motor = motor_rpm(temp,params); % get motor speed
                temp.t = dels_scale(o) / temp.v; % get time taken to complete lap segment
                F_lateral = cornering(temp,curv_scale(o),params); % get the lateral force at the simulated segment of the forward looking simulation
                F_brake = brake(F_lateral,temp,params); %  get the max braking  using the temp future state
                F_aero = aero_force(temp,params);
                F_vehicle = -F_brake - F_aero;
                temp.v = temp.v + temp.t * F_vehicle / params.M; % update temporary velocity using the maximum availible braking force
                temp.F = -F_brake;
                [temp.F_long_load_transfer, temp.a_long] = long_load_transfer(params,tempPrev); % compute bicycle model for this segment
    
                % if the temporary velocity is higher than max velocity at a given
                % track segment (ie going too fast for the next corner) or the braking force availible is small (ie very high lateral loads) apply the brake flag
                % break if this condition is true as there is no need to simulate
                % beyond this point
                if temp.v > velocity_max(o) || F_brake <= 0
                    state.brake_flag = 1;
                    state.brake_index = o;
                    break
                end
                % if the vehicle speed approaches zero remove the brake flag (ie it
                % is too early for the vehicle to start braking) and break
                if temp.v < 1
                    state.brake_flag = 0;
                    break
                end
            end
        end

        % if the driver needs to be braking to make the next corner
        if state.brake_flag == 1
            state.RPM_motor = motor_rpm(state,params); % get motor speed
            state.t = dels_scale(i) / state.v; % get time taken to complete lap segment
            F_lateral = cornering(state,curv_scale(i),params); % get lateral force 
            F_brake = brake(F_lateral,state,params); % get maximum braking force avaible (grip limited)

            F_aero = aero_force(state,params);
            F_vehicle = -F_brake - F_aero;
            state.v = state.v + state.t * F_vehicle / params.M; % compute velocity assuming max braking force applied
            state.F = -F_brake; % braking force made the driving force state
            state.grip_limited = 0; % during braking the vehicle is not grip limited

            % power and regen calculations
            if params.control.regen
                state.T_motor = regen_braking(params, state,F_brake); % get max availible braking torque from motor
            else
                state.T_motor = 0; 
            end
            state.P_motor_drive = state.T_motor * state.RPM_motor * pi / 30; % compute mechanical motor power
            Eff_motor =  motor_efficiency(state.RPM_motor,state.T_motor); % calculate motor efficiency with regen torque
            state.P_motor_draw = state.P_motor_drive / Eff_motor; % compute electrical motor power
            state.P_battery = battery_power(state.P_motor_draw,params,state); % get maximum battery power TODO include inverter efficiency 
            state.I_battery = state.P_battery / pack_voltage(params,state); % compute battery current TODO make pack voltage a funciton of SoC and cell resistance
            state.E = state.P_battery * state.t; % compute pack energy consumed during track segment
        
        else % otherwise accelerating
            state.RPM_motor = motor_rpm(state,params); % get motor speed
            state.v = v_trial; % braking not required so trial velocity is vehicle velocity
            state.F = F_drive; % saved driving force is original driving force
            state.F_veh = F_vehicle;
            state.grip_limited = grip_limited; % save if vehicle is grip limited or not during track segment
            state.T_motor = T_motor; % save motor torque
            state.Eff_motor = Eff_motor; % efficiency of trial drive power is the motor efficincy
            state.P_motor_drive = state.T_motor * state.RPM_motor * pi / 30; % compute motor mechanical power
            state.P_motor_draw = state.P_motor_drive / Eff_motor; % compute motor electrical power 
            state.P_battery = battery_power(state.P_motor_draw, params, state); % get battery power TODO include inverter efficiency 
            state.I_battery = state.P_battery / pack_voltage(params,state); % compute battery current TODO make pack voltage a funciton cell resistance
            state.E = state.P_battery * state.t; % compute pack energy consumed during track segment
            
        end

        % battery model
        state.SoC = update_SoC(params,state);
        state.battery_voltage = pack_voltage(params,state); 
        state.cell_temperature = cell_temperature(params,state); 
        state.cell_losses = cell_losses(params,state);
        storage((lapN-1)*length(curv_scale) + i) = state; % save state structure 
    
        [state.F_long_load_transfer, state.a_long] = long_load_transfer(params,storage); % compute bicycle model for this segment

    end
end

%DATA ANALYSIS
v_data = vertcat(storage.v);
t_lap = sum(vertcat(storage.t)) / Num_Laps
E_endurance = sum(vertcat(storage.E));
E_lap = E_endurance / Num_Laps;
E_endurance_KWh = E_endurance / (3.6 * 10^6)

t_data = cumsum(vertcat(storage.t));

figure;
P_data = 10^-3 * vertcat(storage.P_battery);
P_rms = rms(P_data);
P_max = max(P_data);
plot(P_data)
title("Battery Power Draw Over Lap")
xlabel("Lap Progression")
ylabel("Power from Battery / kW")
annotation('textbox', [0.65, 0.8, 0.1, 0.1], 'String', "RMS Power: "+P_rms+"kW")
annotation('textbox', [0.65, 0.72, 0.1, 0.1], 'String', "Max Power: "+P_max+"kW")
yline(80, 'r:','LineWidth', 2);
ylim([0,120])


I_data = vertcat(storage.I_battery);
I_rms = rms(I_data);

figure;
plot(I_data,"magenta")
title("Battery Current Draw Over Lap")
xlabel("Lap Progression")
ylabel("Current / A")
annotation('textbox', [0.65, 0.8, 0.1, 0.1], 'String', "RMS Current: "+I_rms+"A")

figure;
plot(v_data)
hold on
TF = islocalmin(velocity_max);
plot(x,velocity_max,x(TF),velocity_max(TF),'r*')
xlabel("Lap Progression")
ylabel("Speed / ms^-1")
title("Vehicle Speed vs Max Cornering Speed")
legend("Vehicle Speed","Max Cornering Speed")
ylim([0,50])

figure
yyaxis left
plot(t_data,vertcat(storage.SoC)*100)
ylabel("SoC (%)")
yyaxis right
plot(t_data,vertcat(storage.battery_voltage))
ylabel("Pack OCV (V)")
xlabel("Time (s)")

figure
plot(t_data,vertcat(storage.cell_temperature))
ylabel("Cell Temperature (°C)")
xlabel("Time (s)")

%% thermal sim
T_init = params.ambient_temperature;
heat_cell = (I_data / params.battery.Np).^2 * params.cellR ;
animation_t = [t_data(1:length(curv_scale):end); t_data(end)];

% simulation that computes temp profile at end of each lap using an
% averaged heat flux and average car velocity
[T_animations, model] = thermal_model_2D_transient(T_init,animation_t,mean(heat_cell),mean(v_data));

% simulation that uses the profile to update at the model at every time
% point
% [T_animations, model] = thermal_model_2D_transient_profile(T_init,t_data,heat_cell,v_data,animation_t);


%% Plot at a specific time index (e.g., t = 600 s)
k = length(animation_t)-1;   
figure
pdeplot(model,'XYData',T_animations(:,end),'ColorMap','jet');
xlabel("dimensions (m)")
c = colorbar; c.Label.String = 'Temperature (°C)';
title(sprintf('Temperature at t = %.1f s', t_data(k)));

% Create animated GIF
gifFile = 'thermal_animation.gif';
cbar_interval = 1;
for k = 1:length(animation_t)-1
    pdeplot(model,'XYData',T_animations(:,k),'ColorMap','jet');
    caxis([T_init ceil(max(max(T_animations))/cbar_interval)*cbar_interval])
    c = colorbar; c.Label.String = 'Temperature (°C)';
    title(sprintf('Temperature at t = %.2f s',animation_t(k+1)));
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
