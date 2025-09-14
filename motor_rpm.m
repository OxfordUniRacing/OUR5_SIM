function RPM_motor = motor_rpm(state,params)
    RPM_motor = 60 * state.v * params.gratio / (params.tyre_dia * 25.4 * 10^-3 * pi);
end
