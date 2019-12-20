% 2D FFT method of electromagnetic force and flux density 
% Copyright (c) 2017-2020
% School of Electrical and Electronics of Engineering
% Huazhong University of Science and Technology
% All Rights Reserved
% 
%              Programmer:Lu,Yang yanglu_hust@sina.com 
%              Version:3.0
%              Last modified:11/21/2017
clear;clc;
%% 数据读取与处理
load('Force_36s6p');
p=3;
Qs=36;
Nt=gcd(Qs,p);
n=1500;%speed
f=p*n/60;%fundamental electrical frequency
T=1/f;%period
TimeDomain.Time=linspace(0,T,120);
TimeDomain.Space=linspace(0,350,36);


N.Time=length(TimeDomain.Time);% 时间采样点数
Fs.Time=1/(T/N.Time);% 时间采样频率
N.Space=length(TimeDomain.Space);
Fs.Space=1/(TimeDomain.Space(2)-TimeDomain.Space(1));% 空间采样频率

[Fourier]=fun_fft2d(Frz,N,Fs,f);

View.SpaceOrder=18;
View.TimeOrder=15;
% h = bar3(Fourier.Real.Amplitude(size(Fourier.Real.Amplitude,1)/2+1-View.SpaceOrder:size(Fourier.Real.Amplitude,1)/2+1+View.SpaceOrder,:),0.6);
% stem3(Fourier.Real.Amplitude(:,1:View.TimeOrder),'LineWidth',2);
h=bar3(Fourier.Real.Amplitude(:,1:View.TimeOrder),0.3);
zlim([0 max(max(Fourier.Real.Amplitude))]);
for i = 1:numel(h)
    zData = get(h(i),'ZData');
    zData = repmat(max(zData,[],2),1,4);
    set(h(i),'CData',zData);
    set(h(i),'FaceColor','flat');
end
set(gca,'xlim',[0.5 View.TimeOrder+0.5],'xtick',1:2:View.TimeOrder,'xticklabel',Fourier.Real.TimeOrder(1:2:View.TimeOrder)*f,'Fontsize',12);
set(gca,'ylim',[0.5 View.SpaceOrder+0.5],'ytick',1:2:View.SpaceOrder,'yticklabel',Fourier.Real.SpaceOrder(1:2:View.SpaceOrder),'Fontsize',12);
set(gca,'linewidth',2,'DataAspectRatio', [1 1 10]);
% title('Ridial Force Harmonic','Fontsize',14);
xlabel('Frequency (Hz)','Fontsize',16,'Rotation',20,'FontWeight','bold');
ylabel('Space Harmonic','Fontsize',16,'Rotation',-33,'FontWeight','bold');
zlabel('Force (N)','Fontsize',16,'FontWeight','bold');
grid off
box on
view(3);