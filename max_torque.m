%max_torque.m
%Lewis Blake

function max_torque = max_torque(rpm)
    %peak
    if rpm < 3500
        max_torque = 340;
    end
    if rpm >= 3500
        max_torque = 340 - (rpm-3500)*(340-150)/(8000-3500);
    end
    if rpm > 8000
        max_torque = 0;
    end
    % %continuous
    % if rpm < 6000
    %     max_torque = 50;
    % end
    % if rpm >= 6000
    %     max_torque = 50 - (rpm-6000)*(50-45)/(7000-6000);
    % end
end
