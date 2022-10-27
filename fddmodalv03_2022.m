% FDD modal extraction
% modal damping

% clear all
% clc
% close all
function[fj,svm] = fddmodalv03_2022(nnum,crit)

% nnum =  29 %29 42 % modal data number
%fn = 3.6023%0.354   % natural frequency
Fs = 20;
%nfft=2048*2;
%pcent = 0.5
% load data
sv = load('sv.txt', '-ASCII');
load Uj -mat

fj = sv(:,1);

ncol = length(sv(1,:)');
svj = sv(:,2:ncol);
U1j = Uj(:,1,nnum);
ang = atan2(imag(U1j),real(U1j));%*180/pi;
r1j = abs(U1j);
ms1j = r1j.*cos(ang);

%MAC comparison

nuj = length(Uj);
nms = length(ms1j);
% crit = 0.95;

for j1=1:nuj
    Ujj = Uj(:,:,j1);
    angj = atan2(imag(Ujj),real(Ujj));
    rjj = abs(Ujj);
    msjj = rjj.*cos(angj);
    msjjall(:,:,j1) = msjj;
    for j2=1:nms
        coef=(ms1j'*msjj(:,j2))^2/((ms1j'*ms1j)*(msjj(:,j2)'*msjj(:,j2)));  % MAC
        MAC(j1,j2) = coef;
        if coef>crit
           svm(j1,j2) = svj(j1,j2);
       else
           svm(j1,j2) = 0;
       end
    end
    percent = j1/nuj*100;
    residu = percent/10-floor(percent/10);
    if residu < 0.01
       ['Mode extraction =', num2str(percent),'%']
    end
end

temp = [fj,svm];
save svm02.txt temp -ASCII

% semilogy(fj,svm,'or')
% xlim([0,1])

