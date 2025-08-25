function cell_losses = cell_losses(params,state)
% cell loses due to internal resistance in watts
    cell_losses = (state.I_battery/params.battery.Np)^2 * params.cellR;
end