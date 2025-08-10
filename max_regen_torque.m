%max_regen_torque.m
% Jack Bruford

% Max regen motor torque calculation

function max_regen_torque = max_regen_torque(rpm,params)
    %peak
    power_limit = params.voltage * params.max_charge_Crate * params.pack_Ah;
    pwr_tol = 100;
    max_iterations = 50;
    

    persistent rpm_derate
   
    if(true)
        rpm_min = 0;
        rpm_max = 3100;
    
        for i = 1:max_iterations
            rpm_mid = (rpm_min+rpm_max)/2;
            eff = motor_efficiency(rpm_mid,250);
    
            pwr_est = rpm_mid * 250 * pi/30 / eff;
            bat_pwr = battery_power(pwr_est,params);

            if bat_pwr > power_limit
                rpm_max = rpm_mid; % Update rpm_max for next iteration
            elseif bat_pwr < power_limit - pwr_tol
                rpm_min = rpm_mid; % Update rpm_min for next iteration
            else
                break
            end
        end
        rpm_derate = rpm_mid;
    end

    if rpm > 6100
        max_torque = 0;
    elseif rpm < rpm_derate
        max_torque = 250;
    elseif rpm >= rpm_derate
        % rpm * torque / 9549.3 / eff < power_limit
            
        torque_min = 0;
        torque_max = 250;
        
        for i = 1:max_iterations
            torque_mid = (torque_min + torque_max) / 2;
    
            eff = motor_efficiency(rpm,torque_mid);
            if isnan(eff)
                torque_max = torque_mid;
            else
                pwr_est = rpm * torque_mid * pi/30 / eff;
                bat_pwr = battery_power(pwr_est,params);
                
                
                if bat_pwr > power_limit
                    torque_max = torque_mid;
                elseif bat_pwr < power_limit - pwr_tol
                    torque_min = torque_mid;
                else
                    break
                end
            end
        end

        if isnan(torque_mid)
            keyboard();
        end

        max_torque = torque_mid;
    end

    max_regen_torque = -max_torque;
end
