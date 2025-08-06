%gauss_gen.m
%Lewis Blake

function [x,y]=gauss_gen(r)
x = -(r-1)/2:1:(r-1)/2;
y = normpdf(x,0,r/8);
y = y/((20)*y(ceil(length(y)/2)));

end