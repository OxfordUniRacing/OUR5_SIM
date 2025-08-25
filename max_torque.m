%max_torque.m
%Lewis Blake

% Max motor torque calculation

%Max torque needs to be calculated off of max power with some efficiency
%   losses

% rpm * torque / 9549.3 = power in



function max_torque = max_torque(rpm,params,state)
    %peak
    pwr_tol = 100; %to within how many watts
    max_iterations = 50; %how many iterations for binary search (should never be hit)


    %Check bsaed on battery SOC with max current
    max_current = params.battery.Np * params.cell_Ah * params.max_discharge_Crate;
    max_battery_power = (state.battery_voltage - max_current * params.battery.Ns*params.cellR/params.battery.Np) * max_current;


    %Now derate based on temperature
    % 55deg -> 0.8
    % 58deg -> 0.5
    % 59 deg -> 0.1
    % 60 deg 0 -> 0.01 (keep it moving but should be 0)
    temp_derate = 1.0;

    if params.control.temp_derate
        if state.cell_temperature > 60
            temp_derate = 0.01;
        elseif state.cell_temperature > 59 
            temp_derate = 0.1;
        elseif state.cell_temperature > 58 
            temp_derate = 0.5;
        elseif state.cell_temperature > 55
            temp_derate = 0.8;
        end
    end
    
    power_limit = temp_derate * min(max_battery_power, params.control.max_power);
    


    % See if we need to derate based on voltage
    % Due to internal resistance, we can determine max torque based on
    % current battery voltage and set max battery current 
    % Max battery current should be ~300

    


    rpm_derate = 0;
   
    if(true)
        rpm_min = 0;
        rpm_max = 3100;
    
        for i = 1:max_iterations
            rpm_mid = (rpm_min+rpm_max)/2;
            eff = motor_efficiency(rpm_mid,250);
    
            pwr_est = rpm_mid * 250 * pi/30 / eff;
            %bat_pwr = battery_power(pwr_est,params,state);

            if pwr_est > power_limit
                rpm_max = rpm_mid; % Update rpm_max for next iteration
            elseif pwr_est < power_limit - pwr_tol
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
                %bat_pwr = battery_power(pwr_est,params,state);
                
                
                if pwr_est > power_limit
                    torque_max = torque_mid;
                elseif pwr_est < power_limit - pwr_tol
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
end
