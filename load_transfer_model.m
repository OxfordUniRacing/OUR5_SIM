%load_tranfer_model.m
%Basic point mass sim for determining battery spec
%Lewis Blake

clearvars
close all



addpath("./thermalModel")

%CAR PARAMETERS
params = params();


%TRACK SCALING


%INITIALISATION
state.v = 5;
state.Fz = params.M*params.g;
state.Fz_drive = params.M*params.g*(1-params.M_dist);
state.F = 0;
state.F_long_load_transfer = 0;
state.a_long = 0;
state.SoC = 1.0; % Pack state of charge as a fraction
state.battery_voltage = pack_voltage(params,state);
state.brake_flag = 0;
state.brake_index = 0;
state.cell_temperature = params.ambient_temperature;

%storage(1) = state;

%DERIVED

%SIM

storage = endurance_model(params,state);
%storage = acceleration(params,state);


%% thermal simulation
% run("thermal_simultation.m");
