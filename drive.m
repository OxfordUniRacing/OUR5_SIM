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