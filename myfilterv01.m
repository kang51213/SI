% filter comparison
% Wp = Hz
function [xfilt] = myfilterv01(xt,Fs,Wp,norder)
b = fir1(norder,Wp/(Fs/2),hamming(norder+1));
xfilt = filtfilt(b,1,xt);



