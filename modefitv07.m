% Extracted mode curvefitting
% 07/02/08

% clear all
% clc
% close all

function[w1,d1,rSquare ] = modefitv07(fq,FN,axisNum,fnk)

ftsize = 10;

sdata = load('svm02.txt','-ASCII');
load Uj -mat

f1 = sdata(:,1);
% delf= f1(2); % Hz

col = length(sdata(1,:))-1;  %freq. column excluded
nod = length(sdata(:,1));
svm1 = sdata(:,2:col+1);

% fq = [0.1 0.25];    % 1st [0.07 0.2]; 2nd [1.55, 1.75] frequency range to do curvefitting
% weighting = 1;
% rweig= [0.1, 0.12];  %Hz 1st [1.45, 1.53]; 2nd [1.6, 1.7]
% rweig= [fq(1)*1.3, fq(2)*0.7];  %Hz 1st [1.45, 1.53]; 2nd [1.6, 1.7]

rfq = fq/max(f1);
rndat = ceil(rfq*length(f1));
svmr = svm1(rndat(1):rndat(2),:);
nn = length(svmr(:,1));
fr = f1(rndat(1):rndat(2),:);

%% find nonzero SVmr
% svmrmax = max(svmr')';
kk=0;
for jk=1:length(svmr(:,1))
    tempk = svmr(jk,:);
    nkk = find(tempk ~=0);
    if isempty(nkk) == 0
        kk = kk+1;
        svm(kk,1) = mean(tempk(1,nkk));
        fj(kk,1) = fr(jk,1);
    end
end
% idsv = find(svmrmax ~= 0);  
% idsv = find(svmrMean ~= 0);
% colid = idsv;
% fj = fr(idsv,1);
% % svm = svmrmax(idsv,1);
% svm = svmrMean(idsv,1);
%%


% n1=1
% for j1=1:nn
%     svmr1 = svmr(j1,:);
%     if sum(svmr1) ~=0
%         for j2=1:col
%             if svmr1(j2) ~= 0
%                 fj(n1,1) = fr(j1);
%                 colid(n1,1) = j2;
%                 svm(n1,1) = svmr1(j2);
%                 n1=n1+1;
%                 break
%             end
%         end
%     end
% end
%w1=0.98877;
%d1=0.02;
%m1 = 13820.43107;
%fjj = log(1/m1^2*((fj/w1).^4)./((1-(fj/w1).^2).^2+(2*d1*fj/w1).^2));
logsvm = log(svm);
nd1=length(svm);

% for j1=1:length(svm)
%     if rweig(1)<fj(j1)
%         nii = j1;
%         break
%     end
% end
% 
% for j1=1:length(svm)
%     if rweig(2)<fj(j1)
%         njj = j1;
%         break
%     end
% end


%rweig = rweigf/max(fj);
% nweig = [nii njj];
% 
% weig = ones(nd1,1);
% for jj=nweig(1):nweig(2)
%     weig(jj,1)=weighting;
% end
% figure(1)
% plot(fj,weig)

% f=fittype('log(1/m1^2*((x/w1).^4)./((1-(x/w1).^2).^2+(2*d1*x/w1).^2))'); %'c1*exp(-d1*x)+c2*exp(-d2*x)'
fitFunc = ['log(1/m1^2*((x/',num2str(fnk),').^4)./((1-(x/',num2str(fnk),').^2).^2+(2*d1*x/',num2str(fnk),').^2))']
f=fittype(fitFunc)
options = fitoptions('Method','NonlinearLeastSquares')
options.Normalize = 'off'
options.Robust='LAR'
options.Lower=[0.  0   ]
options.Upper=[0.1 Inf ]
options.StartPoint=[0.02 10]
%options.Algorithm = 'Levenberg-Marquardt'%'Trust-Region'%'Gauss-Newton'
options.Algorithm = 'Trust-Region'
options.DiffMaxChange=1
options.DiffMinChange=1*10^(-24)
options.MaxFunEvals=100000
options.MaxIter=40000
% options.Weights=weig;
[fit1,gof] = fit(fj,logsvm,f,options)

rSquare = gof.rsquare;

figure(13)
plot(fit1,'-r',fj,logsvm)
xlabel('fi(Hz)','FontSize',ftsize,'FontWeight','bold')
ylabel('Selected Singular Value','FontSize',ftsize,'FontWeight','bold')
figname = strcat(FN,'-Fit',num2str(axisNum),'.jpg');
title('Fitting SV Plot','FontSize',ftsize,'FontWeight','bold')
% legend('GHF5R-01-X-tilt','GHF5R-02-X-tilt','GHF5L-01-X-tilt', 'GHF5L-02-X-tilt', 'Location', 'northwest');
set(gca,'FontSize',ftsize,'FontWeight','bold','PlotBoxAspectRatio',[3,1,1])
set(gcf,'position',[300,200,800,400])
print(char(figname),'-djpeg')

%+++++++++++++++++++++++
% draw fredecay
%+++++++++++++++++++++++

%+++++++++++++++++++++++
% draw fredecay
%+++++++++++++++++++++++

d1=fit1.d1;
m1=fit1.m1;
%w1=fit1.w1;
w1=fnk;

% frf = 1/m1^2*((f1/w1).^4)./((1-(f1/w1).^2).^2+(2*d1*f1/w1).^2);
% 
% nfj = length(fj);
% 
% for j1=1:nfj
%     if j1==1
%         ni = 1;
%     else
%         ni = round(fj(j1-1)/delf);
%     end
%     if j1==nfj
%         nj=nod;
%     else
%         nj = round(fj(j1)/delf);
%     end
%     svmj(ni:nj,1:col) = 0;
%     svmj(ni:nj,colid(j1)) = frf(ni:nj,1);
% end
% 
% temp = [f1,svmj];
% 
% save svmj03.txt temp -ASCII
%     
    
    
    
    
    
    
    
    
    
    


