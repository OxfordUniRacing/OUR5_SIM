%motor_efficiency.m

function motor_efficiency = motor_efficiency(rpm,torque)
    torque = abs(torque); % assume symetric motor efficiency map with torque
    persistent F

    if isempty(F)
       S = load('motor_efficiency.mat', 'F');
       F = S.F;
    end

    if rpm > 6100
        motor_efficiency = 0.9;
    else
        motor_efficiency =  F(torque,rpm) / 100;
    end

    if motor_efficiency == 0
        motor_efficiency = 1;
    end
end