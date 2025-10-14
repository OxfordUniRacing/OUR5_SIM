function storage = acceleration_model(params,state)
    time_step = 0.001;
    track_length = 75;


    distance = 0;
    step_number = 1;
    state.t = time_step;

    while distance < track_length 


        state.Fz_drive = params.M * params.g *(1-params.M_dist) + state.F_long_load_transfer;

        [F_drive, T_motor, Eff_motor, grip_limited] = drive(0,state,params);   % returns what the max torque we can get is

        F_aero = aero_force(state,params);
        F_rr = params.M * params.g * params.Crr;

        F_vehicle = F_drive - F_aero - F_rr;
        state.F_areo = F_aero;
        state.F_rr = F_rr;
        
        %Update velocity
        state.v = state.v + time_step * F_vehicle / params.M;
        %Update distance
        distance = distance + state.v * time_step;

        %Update state generally
        state.RPM_motor = motor_rpm(state,params); % get motor speed
        state.v = state.v; 
        state.F = F_drive; % saved driving force is original driving force
        state.F_veh = F_vehicle;
        state.grip_limited = grip_limited; % save if vehicle is grip limited or not during track segment
        state.T_motor = T_motor; % save motor torque
        state.Eff_motor = Eff_motor; % efficiency of trial drive power is the motor efficincy
        state.P_motor_drive = state.T_motor * state.RPM_motor * pi / 30; % compute motor mechanical power
        state.P_motor_draw = state.P_motor_drive / Eff_motor; % compute motor electrical power 
        [state.P_battery,state.I_battery] = battery_power(state.P_motor_draw, params, state); % get battery power TODO include inverter efficiency 
        state.E = state.P_battery * time_step; % compute pack energy consumed during track segment

        storage(step_number) = state;
        step_number = step_number + 1;

        %Update load transfer
        [state.F_long_load_transfer, state.a_long] = long_load_transfer(params,storage);
    end

    time = sum(vertcat(storage.t))
   
    for x = 1:(step_number-1)
        time_data(x) = time_step * (x-1);
    end


    
    figure;
    P_data = 10^-3 * vertcat(storage.P_battery);
    P_rms = rms(P_data);
    P_max = max(P_data);
    plot(time_data,P_data)
    title("Power over time")
    xlabel("Lap Progression")
    ylabel("Power from Battery / kW")
    annotation('textbox', [0.65, 0.8, 0.1, 0.1], 'String', "RMS Power: "+P_rms+"kW")
    annotation('textbox', [0.65, 0.72, 0.1, 0.1], 'String', "Max Power: "+P_max+"kW")
    yline(80, 'r:','LineWidth', 2);
    ylim([0,120])



    v_data = vertcat(storage.v);
    figure;
    plot(time_data,v_data)
    title("Velocity over time")
    

    I_data = vertcat(storage.I_battery);
    I_rms = rms(I_data);
    
    figure;
    plot(time_data,I_data,"magenta")
    title("Battery Current Draw Over Time")
    xlabel("Lap Progression")
    ylabel("Current / A")
    annotation('textbox', [0.65, 0.8, 0.1, 0.1], 'String', "RMS Current: "+I_rms+"A")

end
