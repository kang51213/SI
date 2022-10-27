% force spectrum with pressure analysis for chimney V02
% genralized force spectrum & base shear and moment spectrum
% clear all
% clc
% close all
%function[d1,d2,d3,w1,w2,w3] = mrd_cfv01_alt1(xt,delt)%,   wn,nfft,,pcent,FN,axisNum) ,tp,xtp,rSquare
% smp = 3;
% axis = 'c'
% fname = ['fdecay',axis,'.txt',num2str(smp)];
%%
xt1 = load('xt2.txt','-ASCII');
ndat = floor(length(xt1)/2);
xt_half = xt1(1:ndat,1);
rdat = 1.0;
rn = floor(ndat*rdat);
%%
% t1=fdec(:,1);
nxt = floor(length(xt_half));
t1 = [0:delt*2:(nxt-1)*delt*2]';
tj = t1;
ut1=xt_half;%fdec;
plot(t1,ut1)
%%
% ndat = length(ut1)
for j1=ndat:-1:1
    uti=ut1(j1);
    utj=ut1(j1-1);
    if uti*utj<0
       rdat=j1/ndat;
       break
    end
end

%%

%rn = floor(ndat*rdat);
t = t1(1:rn,1);
ut = ut1(1:rn,1);
%%

f=fittype('c1/sqrt(1-d1^2)*exp(-d1*2*pi*w1*x).*cos(2*pi*sqrt(1-d1^2)*w1*x+z1)+c2/sqrt(1-d2^2)*exp(-d2*2*pi*w2*x).*cos(2*pi*sqrt(1-d2^2)*w2*x+z2)+c3/sqrt(1-d3^2)*exp(-d3*2*pi*w3*x).*cos(2*pi*sqrt(1-d3^2)*w3*x+z3)'); %'c1*exp(-d1*x)+c2*exp(-d2*x)' 
%f=fittype('c1*exp(-d1*2*pi*w1*x).*cos(2*pi*w1*x+z1)+c2*exp(-d2*2*pi*w2*x).*cos(2*pi*w2*x+z2)+c3*exp(-d3*2*pi*w3*x).*cos(2*pi*w3*x+z3)');
options = fitoptions('Method','NonlinearLeastSquares') %NonlinearLeastSquares
options.Normalize = 'off'
options.Robust='LAR'
options.Lower=     [-inf -inf  -inf       0.0     0.0       0     0.1  0.17  0.49    -pi  -pi -pi]
options.Upper=     [ inf  inf   inf       0.05    0.05    0.05    0.15  0.23  0.53     pi   pi  pi]
options.StartPoint=[0.1  0.1  0.1   0.01    0.01    0.01    0.1  0.1  0.1     0    0    0]
%options.Lower=     [-2 -2        0.0    0.0     0.2  0.2      -pi  -pi  ]
%options.Upper=     [ 2  2        0.5    0.5     0.5  0.5       pi   pi  ]
%options.StartPoint=[0.1  0.1     0.06   0.06    0.1  0.1       0    0   ]

%options.Algorithm = 'Levenberg-Marquardt'%'Trust-Region'%'Gauss-Newton'
options.Algorithm = 'Trust-Region'
options.DiffMaxChange=1
options.DiffMinChange=1*10^(-3)
options.MaxFunEvals=100000
options.MaxIter=40000
fit1 = fit(t,ut,f,options)
c1=fit1.c1;
d1=fit1.d1;
w1=fit1.w1;
z1=fit1.z1;

c2=fit1.c2;
d2=fit1.d2;
w2=fit1.w2;
z2=fit1.z2;

c3=fit1.c3;
d3=fit1.d3;
w3=fit1.w3;
z3=fit1.z3;


figure(1)
h=plot(fit1,'-r',t,ut,'ob')
set(h,'MarkerSize',3,'Linewidth',2.1)
    %xlim([min(fj),max(fj)])
    %ylim([0, 25])
    %grid on
    xlabel('Time (s)','FontSize',13,'FontWeight','bold')
    ylabel('Mag.','FontSize',13,'FontWeight','bold')
    figname = ['smp=',smp,'axis=',axis,'drx=',num2str(d1),'dry=',num2str(d2)];
    title(figname,'FontSize',13,'FontWeight','bold')
    %set(gca,'XTick',[0:2:24])
    %set(gca,'YTick',[0:5:25])
    set(gca,'FontSize',11,'FontWeight','bold','PlotBoxAspectRatio',[1.5,1,1])
    set(gcf,'PaperSize',[700,450])
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'position',[300,300,700,450])
    myfile = ['MRD_',axis,num2str(smp),'.tif'];
    print('-dtiff','-r500','-append', myfile)



f11=c1*exp(-d1*2*pi*w1*tj).*cos(2*pi*w1*tj+z1);
f21=c2*exp(-d2*2*pi*w2*tj).*cos(2*pi*w2*tj+z2);
f31=c3*exp(-d3*2*pi*w3*tj).*cos(3*pi*w3*tj+z3);

ffit = f11+f21;


temp = [tj fdec ffit,f11,f21];

figure(2)
plot(tj,ut1,':',tj,f11,'-b',tj,f21,'-k',tj,f31,'-r')

fname=['mrd_results',axis,num2str(smp),'.txt'];

save(fname, 'temp', '-ASCII')
