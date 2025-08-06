%motor_efficiency.m

function motor_efficiency = motor_efficiency(rpm,torque)

if torque >= 90 || torque <= 16
    motor_efficiency = 0.86;
end
if torque > 16 && torque < 90
    if rpm < 1500
        motor_efficiency = 0.86;
    end
    if rpm >= 1100
        motor_efficiency=0.9;
    end
    if rpm >=6000
        motor_efficiency = 0.86;
    end
end

if rpm > 6000 && torque > 16+10*(rpm-6000)/1000
    motor_efficiency = 0.9;
end

if torque >= 22 || torque <= 80
    if rpm >= 1500 || rpm <= 6000
        motor_efficiency = 0.94;
    end
    if rpm >= 6000 && torque > 22+10*(rpm-6000)/1000
        motor_efficiency = 0.94;
    end

    if torque <= 70 || torque >= 26
        if rpm > 2300 && rpm < 6500
            motor_efficiency = 0.95;
        end
    end

    if torque <= 55 && torque >= 32
        if rpm > 2800 && rpm < 5500
            motor_efficiency = 0.96;
        end

    end

end
end