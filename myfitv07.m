% dampin

clear all
clc
close all

adat = load('freedacy.txt','-ASCII');
t1=adat(:,1);
ut1=adat(:,2);
ndat = length(ut1)
rdat = 0.5;
rn = floor(ndat*rdat);
t = t1(1:rn,1);
ut=ut1(1:rn,1);

f=fittype('c1*exp(-d1*w1*x).*cos(w1*x+a1)+c2*exp(-d2*w2*x).*cos(w2*x+a2)'); %'c1*exp(-d1*x)+c2*exp(-d2*x)'
options = fitoptions('Method','NonlinearLeastSquares')
options.Normalize = 'off'
options.Robust='LAR'
%options.Lower=[0 0 0 0 0 0]
%options.Upper=[Inf Inf Inf Inf Inf Inf]
%options.StartPoint=[0.5621 14.12 0.04679 0.09046  0.02073 57.59]
options.Algorithm = 'Levenberg-Marquardt'%'Trust-Region'%'Gauss-Newton'
%options.Algorithm = 'Trust-Region'
options.DiffMaxChange=100000
options.DiffMinChange=1*10^(-24)
options.MaxFunEvals=100000
options.MaxIter=40000
fit1 = fit(t,ut,f,options)
plot(fit1,'-r',t,ut)
