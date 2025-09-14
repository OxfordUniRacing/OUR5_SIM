
function storage = endurance_model(params,state)
    
    load('curve.mat','curv', 'dels','track_length');

    %TRACK PARAMETERS
    endurance_length = 22e3;
    lap_length = endurance_length/params.Num_Laps; %m
    corner_min_rad = 5; %m
    
    curv_scale_val = (1/corner_min_rad) / abs(max(abs(curv))); %min hairpin radius of Xm
    curv_scale = curv * curv_scale_val;
    length_scale = lap_length / track_length;
    dels_scale = dels * length_scale;


    lap_times = zeros(1,params.Num_Laps);
    n_steps = length(curv_scale);

    Fy_max = state.Fz*params.lat_mu;


    % Calculate maximum velocity at every part throughout the track
    x = 1:length(curv_scale);
    
    for k = x
        velocity_max(k) = sqrt(Fy_max/(abs(curv_scale(k)) * params.M));
    end
    
    TF = islocalmin(velocity_max);
    
    %SIM
    
    sim.runover = 50;
    
    for lapN = 1:params.Num_Laps
        for i = 1:n_steps
    
            % make sure that the driver doesnt come off the brake before the
            % corner apex which would cause oscillation
            if (state.brake_flag == 1) && (i < state.brake_index) 
                state.brake_flag = 1;
            else
                state.brake_flag = 0;
            end
    
            state.t = dels_scale(i) / state.v;
            state.Fz_drive = params.M*params.g*(1-params.M_dist) + state.F_long_load_transfer;
            F_lateral = cornering(state,curv_scale(i),params);
            [F_drive, T_motor, efficiency, grip_limited] = drive(F_lateral,state,params);
            F_aero = aero_force(state,params);
    
            T_motor = T_motor * params.control.driver_skill;
    
            F_vehicle = F_drive - F_aero;
            v_trial = state.v + state.t * F_vehicle / params.M;
    
            state.t = dels_scale(i) / state.v; % time taken to travel through the track segment
            state.Fz_drive = params.M * params.g * (1-params.M_dist) + state.F_long_load_transfer; % vertical load on the driven wheels
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
            
                for k = i+1:n_steps+sim.runover
                    % if index exceeds track length loop back to the start of the track
                    if k > n_steps
                        o = k - n_steps;
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
                [state.P_battery, state.I_battery] = battery_power(state.P_motor_draw,params,state); % get maximum battery power TODO include inverter efficiency 
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
                [state.P_battery,state.I_battery] = battery_power(state.P_motor_draw, params, state); % get battery power TODO include inverter efficiency 
                state.E = state.P_battery * state.t; % compute pack energy consumed during track segment
                
            end
    
            % battery model
            state.SoC = update_SoC(params,state);
            state.battery_voltage = pack_voltage(params,state); 
            state.cell_temperature = cell_temperature(params,state); 
            state.cell_losses = cell_losses(params,state);
            storage((lapN-1)*n_steps + i) = state; % save state structure 
        
            [state.F_long_load_transfer, state.a_long] = long_load_transfer(params,storage); % compute bicycle model for this segment
    
        end
        if(state.SoC < 0.1) 
            fprintf("Battery Empty!")
            params.Num_Laps = lapN;
            break;
        end
    end



    v_data = vertcat(storage.v);
    t_lap = sum(vertcat(storage.t)) / params.Num_Laps
    E = sum(vertcat(storage.E));
    E_lap = E / params.Num_Laps;
    E_KWh = E / (3.6 * 10^6);
    E_endurance_KWh = E_KWh
    
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
    xlim([0,947])
    
    figure;
    plot(lap_times)
    hold on
    title("Lap Times over Endurance")
    xlabel("Lap Number")
    ylabel("Lap Time")
    
    
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



end

