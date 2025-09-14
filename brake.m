function F_brake = brake(F_lateral,state,params)
    F_brake = state.Fz * params.long_mu * sqrt(1-(F_lateral/(state.Fz*params.lat_mu)));
    
    if abs(imag(F_brake)) > 0
        F_brake = 0;
    end
    
    F_brake = params.control.driver_skill * F_brake; %  get the max braking force
    F_brake = (1 - params.control.driver_smoothness_alpha) * -state.F + params.control.driver_smoothness_alpha * F_brake;% Smooth braking force (F_brake is positive state.F is negative for braking)

end