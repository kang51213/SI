% resampling of data
function [xtr] = resmpv01(x1,Fso,Fsc)

n1 = Fso/Fsc;
xtr = x1(1:n1:length(x1));



