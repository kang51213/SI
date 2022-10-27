clear all
clc
close all

%% data read
% [FN,PN,FI] = uigetfile('*.txt','MultiSelect', 'on');
% FN = sortrows(FN');
% FN = {'20220906 0340.csv'}
load('x1.mat');
accMerge(:,1) = adata(50000:end,1);
clear adata
load('y1.mat');
accMerge(:,2) = adata(50000:end,1);
clear adata
load('y2.mat');
accMerge(:,3) = adata(50000:end,1);
clear adata
load('z1.mat');
accMerge(:,4) = adata(50000:end,1);
clear adata
%% Input
Fs = 20;
nfft = 8192;%2048*2;
ftsize = 10;
fq = [0.15 0.35
    0.15 0.35];
crit = 0.99;
pcent = 0.8;
%%
for j1=1:1%length(FN)
%     fname = char(FN(j1));
    %fname = '20220905 2350.csv';
%    lmtraw = csvread(fname,1,1);
%     lmtraw = readmatrix(fname,'Range','S:V');
%    araw = lmtraw(:,1:4)*98.1; %gal
    araw = accMerge; %gal
    
    %% FDD
    [temp, aRms(j1,:), aPeak(j1,:)] = fddv04_2022(araw,nfft);
    %% find peaks & natural frequency
    fRng = [0.01,5];
    fi = temp(:,1);
    Sv = temp(:,2:end);
    ind = find(fi>fRng(2));
    Sv(ind,:) = zeros(length(ind),4);
    [pk,lk] = findpeaks(Sv(:,1),'Sortstr','descend');
    fn(j1,:) = sortrows(fi(lk(1:2),1))';
    %% figure
    figure(5)
    semilogy(fi,Sv,'-b',fi(lk(1:2),1),pk(1:2,1),'or');

    xlabel('fi(Hz)','FontSize',ftsize,'FontWeight','bold')
    ylabel('Singular Value','FontSize',ftsize,'FontWeight','bold')
    figname = strcat(FN(j1),'-SV.jpg');
    title('SV Plot','FontSize',ftsize,'FontWeight','bold')
    % legend('GHF5R-01-X-tilt','GHF5R-02-X-tilt','GHF5L-01-X-tilt', 'GHF5L-02-X-tilt', 'Location', 'northwest');
    set(gca,'FontSize',ftsize,'FontWeight','bold','PlotBoxAspectRatio',[3,1,1])
    set(gcf,'position',[300,200,800,400])
    print(char(figname),'-djpeg')
    
    %% save singular value    
    tempj = [fi,Sv];
    save Sv.txt tempj -ASCII
    
    %% modal data
    for j2=1:2
        [fj,svmk] = fddmodalv03_2022(lk(j2),crit);
        [w1Fit(j1,j2),d1Fit(j1,j2),rSquareFit(j1,j2)] = modefitv07(fq(j2,:),FN(j1),j2,fn(j1,j2));
        [w1RD(j1,j2),d1RD(j1,j2),rSquareRD(j1,j2),fn_mrd(j1,j2), damp_mrd(j1,j2),rSquare_mrd(j1,j2)] = dampautov01(fn(j1,j2),Fs,nfft,FN(j1),j2,pcent);
        
        %%
        figure(6)
        semilogy(fi,svmk,'-o');
        xlabel('fi(Hz)','FontSize',ftsize,'FontWeight','bold')
        ylabel('Selected Singular Value','FontSize',ftsize,'FontWeight','bold')
        figname = strcat(FN(j1),'-SV-selected',num2str(j2),'.jpg');
        title('Selected SV Plot','FontSize',ftsize,'FontWeight','bold')
        % legend('GHF5R-01-X-tilt','GHF5R-02-X-tilt','GHF5L-01-X-tilt', 'GHF5L-02-X-tilt', 'Location', 'northwest');
        set(gca,'FontSize',ftsize,'FontWeight','bold','PlotBoxAspectRatio',[3,1,1])
        set(gcf,'position',[300,200,800,400])
        print(char(figname),'-djpeg')

        
    end

    
end

%% Results
dPeak = aPeak./(fn*2*pi()).^2;
dRms = aRms./(fn*2*pi()).^2;
xResults = [dPeak(:,1), dRms(:,1),aPeak(:,1), aRms(:,1), fn(:,1),w1Fit(:,1),w1RD(:,1), fn_mrd(:,1), d1Fit(:,1),d1RD(:,1),damp_mrd(:,1),rSquareFit(:,1), rSquareRD(:,1),rSquare_mrd(:,1)];
yResults = [dPeak(:,2), dRms(:,2),aPeak(:,2), aRms(:,2), fn(:,2),w1Fit(:,2),w1RD(:,2), fn_mrd(:,2), d1Fit(:,2),d1RD(:,2),damp_mrd(:,2),rSquareFit(:,2), rSquareRD(:,2),rSquare_mrd(:,2)];

Varname = {'Dxpeak(cm)' 'Dxrms(cm)' 'Axpeak(gal)' 'Axrms(gal)' 'fnx_si(Hz)' 'fnx_Fit(Hz)' 'fnx_RD(Hz)','fnx_MRD(Hz)','drx_Fit','drx_RD','drx_MRD','R2_Fit','R2_RD','R2_MRD'};
T_SI_x = table(dPeak(:,1), dRms(:,1),aPeak(:,1), aRms(:,1), fn(:,1),w1Fit(:,1),w1RD(:,1), fn_mrd(:,1), d1Fit(:,1),d1RD(:,1),damp_mrd(:,1), rSquareFit(:,1),...
    rSquareRD(:,1),rSquare_mrd(:,1), 'VariableNames',Varname);
Varname = {'Dypeak(cm)' 'Dyrms(cm)' 'Aypeak(gal)' 'Ayrms(gal)' 'fny_si(Hz)' 'fny_Fit(Hz)' 'fny_RD(Hz)','fny_MRD(Hz)','dry_Fit','dry_RD','drx_MRD','R2_Fit','R2_RD','R2_MRD'};
T_SI_y = table(dPeak(:,2), dRms(:,2),aPeak(:,2), aRms(:,2), fn(:,2),w1Fit(:,2),w1RD(:,2), fn_mrd(:,2), d1Fit(:,2),d1RD(:,2),damp_mrd(:,2),rSquareFit(:,2),...
    rSquareRD(:,2),rSquare_mrd(:,2), 'VariableNames',Varname);

writetable(T_SI_x,'acc_results1.xlsx','Sheet','xResults','Range','A1');%,'WriteVariableNames',false);
writetable(T_SI_y,'acc_results1.xlsx','Sheet','yResults','Range','A1');%,'WriteVariableNames',false);


