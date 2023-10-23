function [y] = circmean(x,dim)
%Circular mean; see Batschelet 1981 - JZV
y = angle(sum(exp(1i*x),dim));
end

