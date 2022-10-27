% FDD mode shape generation
% modal damping

clear all
clc
close all

nnum = [29 42]     %263%; % modal data number
%fn = 0.354   % natural frequency
%Fs = 100;
%nfft=2048*2;
%pcent = 0.5
% load data
sv = load('sv.txt', '-ASCII');
load Uj -mat

fj = sv(:,1);

ncol = length(sv(1,:)');
svj = sv(:,2:ncol);
for j1=1:length(nnum)
    U1j = Uj(:,1,nnum(j1));
    ang = atan2(imag(U1j),real(U1j));%*180/pi;
    r1j = abs(U1j);
    ms1j = r1j.*cos(ang);
    msmax = max(abs(ms1j));
    msjj = ms1j/msmax;
    for j2=1:1
        msall(j2,:,j1) = msjj;%msjj(3*(j2-1)+1:3*j2,1);
    end
end

save msall01 msall -mat


