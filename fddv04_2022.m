% FDD

% clear all
% clc
% close all
%%
function[temp,aRms,aPeak] = fddv04_2022(araw,nfft)

%%

%fdir = 'D:\project2007\연구과제\FDDroving\example2\';
%fwname = ['red_k2ascii.txt'];
ncol =4;
Fs = 20;%9.4037; %after resempling
fend = 400; % end of frequency for svd
delt = 1/Fs;
% nfft = 2048*1;
p1 = 20;
nsll = 1;
norder = 50;
Wp = [0.05, 9.9];
Fso = 100;
Fsc = 20;
%adata = load(fwname,'-ASCII');
%%
% lmtraw = csvread('LMT_101F_20220905 2350.csv',0,1);
% araw = lmtraw(:,1:6);
for j1=1:ncol
    araw_filt =  myfilterv01(araw(:,j1),Fs,Wp,norder);
    araw_res(:,j1) = resmpv01(araw_filt,Fso,Fsc);
    
end

adata = araw_res(:,1:ncol);
adata1 = detrend(araw_res(:,1:ncol));
% figure(20)
% plot([1:length(adata)]',adata,'-r',[1:length(adata)]',adata1,'-b');

%% peak, rms acc
aRms = std(adata1);
aPeak = max(abs(adata1));


%%
nsmp = floor(length(adata)/nfft)  % number of sample
dend = nfft*nsmp;
ntime = floor(nsmp/nsll);
xnfj = adata(1:dend,:);   % correction for this example
t = [0:delt:(dend-1)*delt];


for j3=1:ncol
    col = j3;
    xnfi=detrend(xnfj(:,col));   % xnfkk(:,j3) = myfilterv01(xnfi,Fs,Wp,norder);
    xnf1=detrend(xnfj(:,1));
    xnf=xnfi/std(xnf1);
    % PSD through FFT
    window = hanning(nfft);
    noverlap = nfft/2
    dflag='none'
%     [Sx_fft,f]=psd(xnf,nfft,Fs,window,noverlap,dflag);
%     Sx=2*Sx_fft/Fs;
    [P2,f1] = pyulear(xnf,p1,nfft,Fs);
    %PSD2(:,j) = [P2];
    %PSD(:,j) = [Sx];

    PSD2a(:,j3) = [P2];
%     PSDa(:,j3) = [Sx];
    ndd = length(xnf);
   t = [0:delt:delt*(ndd-1)]'; % time
% plot results
    figure(2*j3-1)
    plot(t,xnf)

    xlabel('Time (sec)','FontSize',9,'FontWeight', 'bold')
    ylabel('Acceleration(cm/s^2 )','FontSize',9,'FontWeight', 'bold')
    %title(j3)

    figure(2*j3)
    plot(f1,P2,'r');%,f,Sx,':k');
    xlim([0.01,10]);%max(f1)])
    ylabel('Mag','FontSize',9,'FontWeight', 'bold')
    xlabel('Freq (Hz)','FontSize',9,'FontWeight', 'bold')
    %title(j3)

   
    spl = fix(ndd/nsmp);
    for j4=1:nsmp
       xsmp = xnf((j4-1)*spl+1:j4*spl);
       Amax(j4,j3) = max(xsmp);
       Amin(j4,j3) = min(xsmp);
    end

    %Amax_avg(j3) = mean(Amax);
    %Amin_avg(j3) = mean(Amin);
    A_std(j3) = std(xnf);

end

std1=std(xnf1);  %for normalize



for j5=1:ncol
    for j6=1:ncol
        for j7=1:ntime
            ij = (j7-1)*nfft*nsll+1;
            jj = j7*nfft*nsll;
            xi1=detrend(xnfj(ij:jj,j5)/std1,'constant');
            xj1=detrend(xnfj(ij:jj,j6)/std1,'constant');
            xi=myfilterv01(xi1,Fs,Wp,norder);                 
            xj=myfilterv01(xj1,Fs,Wp,norder);
            [Pi,f]=csd(xi,xj,nfft,Fs,window,noverlap,dflag);
            Pxyi(:,j7) = Pi;
        end
            Pxy(j5,j6,:) = mean(Pxyi')';
    end
    
    percent = j5/ncol*100;
%    if percent > 10 | percent > 20 | percent > 30 | percent > 40 | percent > 50 | percent > 60 | percent > 70 | percent > 80 | percent ==100
      ['csd calculation =', num2str(percent),'%']
      %    end

    
end


for j=1:length(P2)
    temp = Pxy(:,:,j);
    [U,S,V] = svd(temp);
    Sj(:,j) = diag(S);
    Uj(:,:,j) = U;   % singular vectors
    if f(j) > fend
        break
    end
    
    percent = j/length(P2)*100;
    %if percent > 10 | percent > 20 | percent > 30 | percent > 40 | percent > 50 | percent > 60 | percent > 70 | percent > 80 | percent ==100
    ['singular calculation =', num2str(percent),'%']
    %end
    
end

Sv = Sj'; % singular values

'save Sv and Uj'

temp = [f(1:j),Sv];
% save Sv.txt temp -ASCII
save Uj Uj -mat;

%%
% figure
% semilogy(f,Sv)
% xlim([0,2])

