function params = params()

    params.M = 320; %kg
    params.M_dist = 0.40; %distribution mass over front wheels
    params.gratio = 4.0; %gear reduction ratio

    params.COG_h = 0.3; %m
    params.wheelbase = 1.525; %m

    %Tyres
    params.lat_mu = 1.2; %lateral tyre coeff. friction
    params.long_mu = 1.4; %longitudinal tyre coeff. friction
    params.tyre_dia = 16; %tyre diameter, inches
    params.Crr = 0.012;
    
    %aero
    params.Cd= 1; % drag co-efficient
    params.Cl = 0.1;% lift co-efficient (NOT IMPLEMENTED)
    params.frontal_area = 1.2; %m^2
    params.air_density = 1.225; 
    
    % BATTERY PARAMETERS
    params.max_charge_Crate = 1; % max charging C-rate
    params.max_discharge_Crate = 10; % max charging C-rate 
    params.battery.Np = 6; % number of cells in parallel in a cell group
    params.cell_Ah = 4.5; % battery pack Amp-hours
    params.pack_Ah = params.cell_Ah * params.battery.Np; % battery pack Amp-hours
    params.cellR = 15e-3; % cell resistance 
    params.cellV = 4.2;
    params.battery.Ns = 90;
    params.battery.cell_specific_heat = 830; % 830 typical for a NCA cell, 1040 would be typical for NMC
    params.battery.cell_mass = 70e-3; % from molicell datasheet, 70g
    params.battery.cellRth = 61;
    %EFFICIENCIES
    %Motor efficiency is in seperate "motor_efficiecy.m"
    params.efficiency.mechanical = 0.92;
    
    %CONTROL
    params.control.driver_skill = 1.0; %Driver skill factor (~0.5 to 1), acts as derate
    params.control.driver_smoothness_alpha = 0.95; % smoothing factor, 0 = slow change, 1 = instant change
    params.control.max_power = 60e3;
    params.control.regen = true; % toggle regen on or off
    params.control.temp_derate = false;
    
    %ENVIRONMENT
    
    params.ambient_temperature = 25;
    params.g = 9.81;


    %ENDURANCE TRACK
    params.Num_Laps = 27;

end