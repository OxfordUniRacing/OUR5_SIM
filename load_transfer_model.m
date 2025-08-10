%load_tranfer_model.m
%Basic point mass sim for determining battery spec
%Lewis Blake

clearvars
close all

load('curve.mat','curv', 'dels','track_length');

%CAR PARAMETERS
params.M = 320; %kg
params.M_dist = 0.5; %dist mass over front wheels
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
params.max_charge_Crate = 3; % max charging C-rate
params.max_discharge_Crate = 10; % max charging C-rate 
params.battery.Np = 5; % number of cells in parallel in a cell group
params.cell_Ah = 4.5; % battery pack Amp-hours
params.pack_Ah = params.cell_Ah * params.battery.Np; % battery pack Amp-hours
params.cellR = 15e-3; % cell resistance 
params.cellV = 2.4;
params.battery.Ns = 90;

%EFFICIENCIES
%Motor efficiency is in seperate "motor_efficiecy.m"
params.efficiency.mechanical = 0.95;

%CONTROL
params.control.driver_skill = 0.95; %Driver skill factor (~0.5 to 1), acts as derate
params.control.regen = true; % toggle regen on or off

%ENVIRONMENT
g = 9.81;

%TRACK PARAMETERS
lap_length = 1000; %m
corner_min_rad = 5; %m
Num_Laps = 22;

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
        T_motor = F_drive * (params.tyre_dia * 25.4 * 10^-3 / 2) / (params.gratio * params.efficiency.mechanical);
        grip_limited = 1;
    else
        F_drive = F_wheel_max;
        T_motor = T_motor_max;
        grip_limited = 0;
    end
    efficiency = motor_efficiency(RPM_motor,T_motor);
end

function [F_long_load_transfer, a_long] = long_load_transfer(params,storage)
    a_long = storage(end).F / params.M;
    F_long_load_transfer = params.M * a_long * params.COG_h / params.wheelbase;
end

%SIM
for lapN = 1:Num_Laps
    for i = 1:length(curv_scale)
    
        state.brake_flag = 0;
    
        state.t = dels_scale(i) / state.v;
        state.Fz_drive = params.M*g*(1-params.M_dist) + state.F_long_load_transfer;
        F_lateral = cornering(state,curv_scale(i),params);
        [F_drive, T_motor, efficiency, grip_limited] = drive(F_lateral,state,params);
        aero_force = 0.5 * params.air_density * (state.v^2) * params.Cd * params.frontal_area;
        F_drive = F_drive * params.control.driver_skill - aero_force;
        T_motor = T_motor * params.control.driver_skill;
        temp.v = state.v + state.t * F_drive / params.M;
        v_trial = temp.v;
    
        state.brake_flag = 0; 
    
    
        state.t = dels_scale(i) / state.v; % time taken to travel through the track segment
        state.Fz_drive = params.M * g * (1-params.M_dist) + state.F_long_load_transfer; % vertical load on the driven wheels
        F_lateral = cornering(state, curv_scale(i), params); % lateral load on the tyres
        [F_drive, T_motor, Eff_motor, grip_limited] = drive(F_lateral,state,params); % drive function returns key variables for ideal drive 
        F_drive = F_drive * params.control.driver_skill; % scale by driver skill fraction parameter
        T_motor = T_motor * params.control.driver_skill; % scale by driver skill fraction parameter
        
        temp.v = state.v + state.t * F_drive / params.M; % temporary vehicle velocity
        
        v_trial = temp.v; 
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
            
            t = dels_scale(o) / temp.v; % get time to complete the next track segment 
            F_lateral = cornering(temp,curv_scale(o),params); % get the lateral force at the simulated segment of the forward looking simulation
            F_brake = params.control.driver_skill * brake(F_lateral,state,params); %  get the max braking force
            temp.v = temp.v - t * F_brake / params.M; % update temporary velocity using the maximum availible braking force
            % if the temporary velocity is higher than max velocity at a given
            % track segment (ie going too fast for the next corner) or the braking force availible is small (ie very high lateral loads) apply the brake flag
            % break if this condition is true as there is no need to simulate
            % beyond this point
            if temp.v > velocity_max(o) || abs(F_brake) < 1000
                state.brake_flag = 1;
                break
            end
            % if the vehicle speed approaches zero remove the brake flag (ie it
            % is too early for the vehicle to start braking) and break
            if temp.v < 1
                state.brake_flag = 0;
                break
            end
        end
        
        % if the driver needs to start braking to make the next corner
        if state.brake_flag == 1
            state.RPM_motor = motor_rpm(state,params); % get motor speed
            state.t = dels_scale(i) / state.v; % get time taken to complete lap segment
            F_lateral = cornering(state,curv_scale(i),params); % get lateral force 
            F_brake = params.control.driver_skill * brake(F_lateral,state,params); % get maximum braking force avaible (grip limited)
            state.v = state.v - state.t * F_brake / params.M; % compute velocity assuming max braking force applied
            state.F = -F_brake; % force on vehicle is equal to the braking force
            state.grip_limited = 0; % during braking the vehicle is not grip limited
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
        else
            state.RPM_motor = motor_rpm(state,params); % get motor speed
            state.v = v_trial; % braking not required so trial velocity is vehicle velocity
            state.F = F_drive; % set vehicle force to driving force
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
    
        storage((lapN-1)*length(curv_scale) + i) = state; % save state structure 
    
        [state.F_long_load_transfer, state.a_long] = long_load_transfer(params,storage); % compute bicycle model for this segment
    
    end
end

%DATA ANALYSIS
v_data = vertcat(storage.v);
t_lap = sum(vertcat(storage.t)) / Num_Laps;
E_endurance = sum(vertcat(storage.E));
E_lap = E_endurance / Num_Laps;
E_endurance_KWh = E_endurance / (3.6 * 10^6)

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

 