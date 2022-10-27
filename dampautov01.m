% damping calculation with auto spectrum
% clear all
% clc
% close all
function[fn,damp,rSquare,fn_mrd, damp_mrd,rSquare_mrd] = dampautov01(fn,Fs,nfft,FN,axisNum,pcent)

%nnum = 307%263%; % modal data number
% fn = 0.1367   % natural frequency
% Fs = 20;
% nfft=2048*2; % <=== don't forget

% load data
svm1 = load('svm02.txt', '-ASCII');% <=== don't forget'svm1_modify_8192_05.txt'
nkk = length(svm1(1,:));
svm = svm1(:,2:nkk);
load Uj -mat

% inverse fft
nuj = length(Uj);

for j=1:nuj
    Vi = Uj(:,:,j);
    crit=sum(svm(j,:));
    if crit==0
        svr(:,:,j) =zeros(nkk-1,nkk-1);
    else
        svr(:,:,j) = Vi*diag(svm(j,:))*Vi';
    end
    
    percent = j/nuj*100;
    residu = percent/10-floor(percent/10);
    if residu < 0.05
       ['Auto spectrum generation =', num2str(percent),'%']
    end
    
end

nsv = length(svm(:,1));
for j6=1:nsv
    Sjji(:,j6) = diag(svr(:,:,j6));
end

Sjj(:,1) = sum(Sjji)';%svr(1,1,:)+svr(2,2,:);

xti = ifft(Sjj);
ang1 = atan2(imag(xti),real(xti));
xt = abs(xti).*cos(ang1);

% save xt2.txt xt -ASCII

% damping evaluation
'Damping evaluation~'
% pcent = 0.8
[damp,tp,xtp,rSquare] = polydampv01(xt,fn*2*pi,nfft,1/Fs,pcent,FN,axisNum);
[d1,d2,w1,w2,rSquare_mrd] = mrd_cfv01_alt1_2022(xt,1/Fs,FN,axisNum,pcent);
if axisNum ==1
    damp_mrd = d1;
    fn_mrd = w1;
else
    damp_mrd = d2;
    fn_mrd = w2
end
    
nxt = length(xt);
delt = 1/Fs;
t = [0:delt*2:(nxt-1)*delt*2]';

temp =[t(1:floor(nxt/2)),xt(1:floor(nxt/2),1)];
save freedacy.txt temp -ASCII
clear temp 
temp = [tp,xtp];
save peak.txt temp -ASCII
