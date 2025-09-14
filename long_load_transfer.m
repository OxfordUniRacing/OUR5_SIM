function [F_long_load_transfer, a_long] = long_load_transfer(params,storage)
    a_long = storage(end).F / params.M;
    F_long_load_transfer = params.M * a_long * params.COG_h / params.wheelbase;
end