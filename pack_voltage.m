function pack_voltage = pack_voltage(params,state)
% fit to cell data in reference folder
% Polynomial Curve Fit (poly7)
% f(x) = p1*x^7 + p2*x^6 + p3*x^5 + p4*x^4 + p5*x^3 + p6*x^2 + p7*x + p8
% 
% Coefficients and 95% Confidence Bounds
%        Value      Lower      Upper  
% p1    -0.0051    -0.0062    -0.0040
% p2    0.0781     0.0598     0.0964 
% p3    -0.4777    -0.5976    -0.3578
% p4    1.4651     1.0665     1.8637 
% p5    -2.3326    -3.0394    -1.6257
% p6    1.7759     1.1340     2.4178 
% p7    -0.7824    -1.0413    -0.5235
% p8    5.0095     4.9763     5.0428 
% 
% Goodness of Fit
%              Value  
% SSE         0.0271 
% R-square    0.9993 
% DFE         78.0000
% Adj R-sq    0.9992 
% RMSE        0.0186 
p = [-0.005075, 0.07809, -0.4777, 1.465, -2.333, 1.776, -0.7824, 5.01];

AhRemaining = state.SoC * params.cell_Ah;

cellV = polyval(p, params.cell_Ah - AhRemaining);

pack_voltage = cellV * params.battery.Ns;
