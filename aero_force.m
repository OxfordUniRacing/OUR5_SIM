function F_aero = aero_force(state,params)
    F_aero = 0.5 * params.air_density * (state.v^2) * params.Cd * params.frontal_area;
end