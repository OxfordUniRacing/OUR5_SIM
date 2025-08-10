function [T_motor, F_regen] = regen_braking(params,state,F_brake)
    % Function to calculate the regen torque for a given state and braking
    % required braking force. Regen torque will likely not provice all required
    % braking force it is assumed the remaining torque is supplied by the
    % mechanical brakes

    T_wheel = -F_brake * params.tyre_dia * 25.4 * 10^-3 / 2; % torque at the wheel
    T_motor = T_wheel / (params.gratio * params.efficiency.mechanical); % torque at the motor
    
    T_motor = max(T_motor, max_regen_torque(state.RPM_motor,params,state)); % limits the negative torque to the max that can be supplied or is required

end