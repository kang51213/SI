%%
clear all
clc
close all

%% data read
[FN,PN,FI] = uigetfile('*.csv','MultiSelect', 'on');
FN = sortrows(FN');

%%
ncol =2;
Fs = 100;%9.4037; %after resempling
fend = 400; % end of frequency for svd
delt = 1/Fs;
nfft = 2048*2;
p1 = 20;
nsll = 1;
norder = 50;
Wp = [0.05, 9.9];
Fso = 100;
Fsc = 100;
%adata = load(fwname,'-ASCII');
%% accelearation
for j1=1:length(FN)
    fname = char(FN(j1)); 
    lmtraw = csvread(fname,0,1);
    araw = lmtraw(:,1:6)*98.1; %gal
    for j2=1:6
        araw_filt =  detrend(myfilterv01(araw(:,j2),Fs,Wp,norder));
        araw_res(:,j2) = resmpv01(araw_filt,Fso,Fsc);
     %  araw_res(:,j2) = detrend(resmpv01(araw(:,j2),Fso,Fsc));
    end
    adata = araw_res(:,1:2);
    vdata = lmtraw(:,7)*20;
    tj = [0:delt:(length(adata)-1)*delt]';
    figure(1)
    plot(tj,adata(:,1),'-r',tj,adata(:,2),'-b')
    figure(2)
    plot(tj,vdata)
    peak_adata(j1,1:2) = max(abs(adata));
    vdata10(j1,1) = max(vdata);
    clear araw_filt araw_res adata
end
%%
figure(3)
tjj = [0:10:10*(length(peak_adata(:,1))-1)]';
plot(tjj,peak_adata(:,1),'or')
ylim([0,5])

figure(4)
plot(tjj,peak_adata(:,2),'^b')
ylim([0,5])

figure(5)
plot(tjj,vdata10,'ob')
%ylim([0,10])

figure(6)
plot(vdata10, peak_adata(:,1),'ob')
ylim([0,5])
xlim([0,20])

figure(7)
plot(vdata10, peak_adata(:,2),'ob')
ylim([0,5])
xlim([0,20])

%% velocity


