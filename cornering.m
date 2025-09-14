function F_lateral = cornering(state,curv,params)
    F_lateral = params.M * state.v^2 * abs(curv);
end