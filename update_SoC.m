function SoC = update_SoC(params,state)
    
    SoC = state.SoC - state.I_battery * state.t / params.pack_Ah / 3600;
end
    
