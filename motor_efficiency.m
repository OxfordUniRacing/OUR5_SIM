%motor_efficiency.m

% efficiency of motor + inverter

function motor_efficiency = motor_efficiency(rpm,torque)
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