function cell_temperature = cell_temperature(params,state)
    cell_thermal_mass = params.battery.cell_specific_heat * params.battery.cell_mass;
    
    Q_cell = (state.I_battery/params.battery.Np)^2 * params.cellR;
    Q_cell_cooling = 0; % initial assumption cell is insulated
    cell_temperature = state.cell_temperature + Q_cell/cell_thermal_mass;

end