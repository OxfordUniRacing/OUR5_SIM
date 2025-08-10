function battery_power = battery_power(motor_power,params,state)
    % If there is a certain power request
    % Assume pack voltage of 90*3.6V = 324V
    % cells have internal resistance of 0.013 ohms
    % current split 5 ways
    % (V_pack - V_drop) * Current = motor_power
    % (V_pack - R_cell * current/5) * Current = motor_power
    % quadratic relation of form ax^2 + bx + c = 0
    % a = R_cell/25
    % b = V_pack
    % c = - motor_power

    a = -params.cellR / params.battery.Np; % cell group resistance
    b = pack_voltage(params,state);
    c = - motor_power;

    if isnan(c)
        keyboard();
    end


    current = roots([a b c]);

    min_current = -params.max_charge_Crate * params.pack_Ah;
    max_current = params.max_discharge_Crate * params.pack_Ah;

    mask = current > min_current & current < max_current;

    battery_power = current(mask) * state.battery_voltage;
end