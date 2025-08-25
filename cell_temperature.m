function cell_temperature = cell_temperature(params,state)
    cell_thermal_mass = params.battery.cell_specific_heat * params.battery.cell_mass;
    
    Q_cell = cell_losses(params,state); % cell losses (W)
    Q_cell_cooling = (state.cell_temperature - params.ambient_temperature) / params.battery.cellRth; % cell cooling (W)
    cell_temperature = state.cell_temperature + (Q_cell - Q_cell_cooling)/cell_thermal_mass;

end