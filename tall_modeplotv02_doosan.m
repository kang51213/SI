% mode plot for tall building
% 4->y1, 5->y2, 6->x1


clear all
clc
close all

B=30;
D=30;
fwname = 'msall'
nset=1;
numsensor=1
num = 2 %<== 모드번호
mag = 10
story = [1];
H=[40 ];


n=1;
for j1=1:nset
    if j1<10
        id = ['0',num2str(j1)];
    else
        id = num2str(j1);
    end
    fname = [fwname,id];
    load(fname,'-mat');
    
    msall = [msall(:,1,:),msall(:,2,:),msall(:,2,:)];  
  %  msall = [msall(:,1,:),msall(:,2,:),msall(:,4,:)];  
    
    temp(:,:) = msall(:,:,num);
    if j1==1
        ms1 = [story(n,:)' temp];
    else
        ms11 = [story(n,:)' temp];
        ms1 = [ms1;ms11(2:numsensor,:)];
    end
    clear temp
    n=n+1;
    
end

r = length(ms1(:,1));
c = length(ms1(1,:));
ms2 = zeros(r,c);
clear temp
for j1=1:r
    maxfl = max(ms1(:,1))
    for j2=1:r
        if ms1(j2,1)==maxfl
            ms2(j1,:)=ms1(j2,:);
            ms1(j2,:)=zeros(1,c);
        end
    end
end

for j1=1:r
    ms3(r-(j1-1),:)=ms2(j1,:);
end

ms = [ms3(:,1),H,ms3(:,2:4)];%load('storymode4.txt','-ASCII');

H = ms(:,2);
nstory = length(ms(:,1));
x = ms(:,3)*mag;
y1 = ms(:,4)*mag;
y2 = ms(:,5)*mag;

c1(1,:) = [0,0,0];
c2(1,:) = [B,0,0];
c3(1,:) = [B,D,0];
c4(1,:) = [0,D,0];

d1(1,:) = [0,0,0];
d2(1,:) = [B,0,0];
d3(1,:) = [B,D,0];
d4(1,:) = [0,D,0];

cr = [c1;c2;c3;c4;c1];

cm = cr;
for j1=1:nstory
    del = y1(j1)-y2(j1);
    r1 = [0,0,H(j1)]+[x(j1),y1(j1),0];
    r2 = [B,0,H(j1)]+[x(j1),y2(j1),0];
    r3 = r2+[del*D/B,sqrt(-del^2+D^2),0];
    r4 = r1+[del*D/B,sqrt(-del^2+D^2),0];
    m1 = [0,0,H(j1)];
    m2 = [B,0,H(j1)];
    m3 = [B,D,H(j1)];
    m4 = [0,D,H(j1)];
    
       cr = [cr;[r1;r2;r3;r4;r1]];
       cm = [cm;[m1;m2;m3;m4;m1]];
       
   c1(j1+1,:) = r1;
   c2(j1+1,:) = r2;
   c3(j1+1,:) = r3;
   c4(j1+1,:) = r4;
   d1(j1+1,:) = m1;
   d2(j1+1,:) = m2;
   d3(j1+1,:) = m3;
   d4(j1+1,:) = m4;

end

h=plot3(cr(:,1),cr(:,2),cr(:,3),'-ob',c2(:,1),c2(:,2),c2(:,3),'-b',c3(:,1),c3(:,2),c3(:,3),'-b',c4(:,1),c4(:,2),c4(:,3),'-b',cm(:,1),cm(:,2),cm(:,3),':r',d2(:,1),d2(:,2),d2(:,3),':r',d3(:,1),d3(:,2),d3(:,3),':r',d4(:,1),d4(:,2),d4(:,3),':r');
axis([-50 50   -50 50  0 80])
set(gca,'FontSize',9,'FontWeight','bold','PlotBoxAspectRatio',[1,1,3])
set(gcf,'position',[100,100,500,800])
set(h,'MarkerSize',7,'Linewidth',2.1)




