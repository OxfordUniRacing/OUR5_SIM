function pack_voltage = pack_voltage(params,state)
% fit to cell data in reference folder
% Fit Name: untitled fit 1
% 
% Polynomial Curve Fit (poly7)
% f(x) = p1*x^7 + p2*x^6 + p3*x^5 + p4*x^4 + p5*x^3 + p6*x^2 + p7*x + p8
% 
% Coefficients and 95% Confidence Bounds
%        Value      Lower      Upper  
% p1    -0.0071    -0.0086    -0.0055
% p2    0.0980     0.0750     0.1209 
% p3    -0.5393    -0.6747    -0.4040
% p4    1.4887     1.0837     1.8937 
% p5    -2.1331    -2.7795    -1.4867
% p6    1.4616     0.9333     1.9899 
% p7    -0.5795    -0.7713    -0.3878
% p8    4.1730     4.1509     4.1952 
% 
% Goodness of Fit
%              Value  
% SSE         0.0120 
% R-square    0.9993 
% DFE         78.0000
% Adj R-sq    0.9992 
% RMSE        0.0124 
p = [-0.0071, 0.0980, -0.5393, 1.4887, -2.1331, 1.4616, -0.5795, 4.1730];

AhRemaining = state.SoC * params.cell_Ah;

cellV = polyval(p, params.cell_Ah - AhRemaining);

pack_voltage = cellV * params.battery.Ns;
