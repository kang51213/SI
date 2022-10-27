%damping evaluation using polyfit
%1-dof 
% pcent: range of damping evaluation from free decay

function[damp,tp,xtp,rSquare] = polydampv01(xt,wn,nfft,delt,pcent,FN,axisNum)

ftsize = 10;

% find (+) peak 
nxt = floor(length(xt)/2*pcent);
t = [0:delt*2:(nxt-1)*delt*2]';
xth = xt(1:nxt,1);

n1=1

for j2=2:nxt-1
    xt1 = xth(j2)-xth(j2-1);
    xt2 = xth(j2+1)-xth(j2);
    temp = xt1*xt2;
    
    if xt1>0 & temp<0
        tp(n1,1) = t(j2,1);
        xtp(n1,1) = xth(j2,1);
        n1=n1+1;
    end
end

yi = log(xtp);
[a,S] = polyfit(tp,yi,1);
rSquare = 1 - S.normr^2 / norm(yi-mean(yi))^2;
% P = polyfitn(tp,yi,1);
% rSquare = P.R2;
% a = P.p;
yfit = a(1)*tp+a(2);
figure(11)
plot(tp,yi,'o',tp,yfit)
yexp = exp(a(2))*exp(a(1)*t);
figure(12)
plot(t,xth,t,yexp,'-r')
xlabel('Time(s)','FontSize',ftsize,'FontWeight','bold')
ylabel('Amplitude','FontSize',ftsize,'FontWeight','bold')
figname = strcat(FN,'-Free Decay',num2str(axisNum),'.jpg');
title('Free Decay Plot','FontSize',ftsize,'FontWeight','bold')
% legend('GHF5R-01-X-tilt','GHF5R-02-X-tilt','GHF5L-01-X-tilt', 'GHF5L-02-X-tilt', 'Location', 'northwest');
set(gca,'FontSize',ftsize,'FontWeight','bold','PlotBoxAspectRatio',[3,1,1])
set(gcf,'position',[300,200,800,400])
print(char(figname),'-djpeg')



damp = -a(1)/wn



