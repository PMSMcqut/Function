function [TimeDomain,Row,Column,Fourier,View]=ForceFluxCalculation(TimeDomain,N,Fs)
%% 电磁力计算
u0=4*pi*10^-7;
% TimeDomain.Fr=(TimeDomain.Br.^2-TimeDomain.Bt.^2)./(2*u0);
TimeDomain.Fr=(TimeDomain.Br.^2)./(2*u0);
% TimeDomain.Fr=(-TimeDomain.Bt.^2)./(2*u0);
% TimeDomain.Ft=TimeDomain.Br.*TimeDomain.Bt./u0;
[Row.F,Column.F]=size(TimeDomain.Fr);
%% 电磁力信号FFT分析
% 径向电磁力分解
% Fourier.Fr.ComplexData=fft2(TimeDomain.Fr,N.Space,N.Time);
Fourier.Fr.ComplexData=fft2(TimeDomain.Fr);
Fourier.Fr.Amplitude=4*abs(Fourier.Fr.ComplexData)/(N.Space*N.Time);% 幅值转化
Fourier.Fr.Amplitude(1,:)=Fourier.Fr.Amplitude(1,:)/2;% 直流分量幅值转化
Fourier.Fr.Amplitude(:,1)=Fourier.Fr.Amplitude(:,1)/2;% 直流分量幅值转化
Fourier.Fr.Amplitude=fftshift(Fourier.Fr.Amplitude);
Fourier.Fr.ComplexData=fftshift(Fourier.Fr.ComplexData);
Fourier.Fr.Phase=angle(Fourier.Fr.ComplexData)*180/pi;% 相位求解
% 切向电磁力分解
% Fourier.Ft.ComplexData=fft2(TimeDomain.Ft,N.Space,N.Time);
% Fourier.Ft.Amplitude=4*abs(Fourier.Ft.ComplexData)/(N.Space*N.Time);% 幅值转化
% Fourier.Ft.Amplitude(1,:)=Fourier.Ft.Amplitude(1,:)/2;% 直流分量幅值转化
% Fourier.Ft.Amplitude(:,1)=Fourier.Ft.Amplitude(:,1)/2;% 直流分量幅值转化
% Fourier.Ft.Amplitude=fftshift(Fourier.Ft.Amplitude);
% Fourier.Ft.ComplexData=fftshift(Fourier.Ft.ComplexData);
% Fourier.Ft.Phase=angle(Fourier.Ft.ComplexData)*180/pi;% 相位求解

Fourier.Force.Frequency=(0:1:N.Time/2)*Fs.Time/N.Time;% 时间频率转化
Fourier.Force.TimeOrder=Fourier.Force.Frequency/Fourier.Force.Frequency(2);% 时间阶次求解
Fourier.Force.SpaceOrder=(0:1:N.Space/2)*Fs.Space/N.Space;
Fourier.SpaceOrder=Fourier.Force.SpaceOrder/Fourier.Force.SpaceOrder(2);% 空间阶次求解

%% 图形绘制
% 电磁力分析结果显示
set(0,'defaultfigurecolor','w')
figure(1)
View.Force.SpaceOrder=48;
View.Force.TimeOrder=60;

% 径向电磁力时域图
s(1)=subplot(2,2,1);
surface(TimeDomain.Time,TimeDomain.ElecAngle,TimeDomain.Fr);
shading interp;
axis([0 0.008 0 360 0 500000]);
title(s(1),'Radial force density','Fontsize',14);
zlabel('\itP\rm/(N/m^2)','Fontsize',14);
set(gca,'xtick',0:0.002:0.008);
set(gca,'xticklabel',0:2:8);
set(gca,'ytick',0:40:360);
set(gca,'yticklabel',0:40:360);
view(3);
box on
grid on
% 
% % 切向电磁力时域图
% s(2)=subplot(2,2,2);
% surface(TimeDomain.Time,TimeDomain.ElecAngle,TimeDomain.Ft);
% shading interp;
% axis([0 0.008 0 360 0 500000]);
% title(s(2),'Tangential force density','Fontsize',14);
% zlabel('\itP\rm/(N/m^2)','Fontsize',14);
% set(gca,'xtick',0:0.002:0.008);
% set(gca,'xticklabel',0:2:8);
% set(gca,'ytick',0:40:360);
% set(gca,'yticklabel',0:40:360);
% view(3);
% box on
% grid on
% 
% 径向电磁力幅频图
s(3)=subplot(2,2,3);
View.Fr.Amplitude=Fourier.Fr.Amplitude(floor(Row.F/2)+1:floor(Row.F/2)+1+View.Force.SpaceOrder,floor(Column.F/2)+1-View.Force.TimeOrder/2:floor(Column.F/2)+1+View.Force.TimeOrder/2);
%View.Fr.Amplitude=Fourier.Fr.Amplitude(floor(Row.F/2)+1-4:floor(Row.F/2)+1+4,floor(Column.F/2)+2:floor(Column.F/2)+1+View.Force.TimeOrder/2);
h = bar3(View.Fr.Amplitude,0.6);
zlim([0 200000]);
for i = 1:numel(h)
    zData = get(h(i),'ZData');
    zData = repmat(max(zData,[],2),1,4);
    set(h(i),'CData',zData);
    set(h(i),'FaceColor','flat');
end
%imagesc(View.Fr.Amplitude);
set(gca,'xtick',1:5:View.Force.TimeOrder+1);
set(gca,'xticklabel',-View.Force.TimeOrder/2:5:View.Force.TimeOrder/2);
title(s(3),'Radial Force Density fft','Fontsize',14);
% xlabel('Time Harmonic','Fontsize',14);
set(gca,'ytick',1:5:View.Force.SpaceOrder+1);
set(gca,'yticklabel',0:5:View.Force.SpaceOrder);
% ylabel('Space Harmonic','Fontsize',14);
zlabel('\itP\rm/(N/m^2)','Fontsize',14);
set(gca,'linewidth',1);
grid on
box on
view(3);
% 
% % 切向电磁力幅频图
% s(4)=subplot(2,2,4);
% View.Ft.Amplitude=Fourier.Ft.Amplitude(floor(Row.F/2)+1:floor(Row.F/2)+1+View.Force.SpaceOrder,floor(Column.F/2)+1-View.Force.TimeOrder/2:floor(Column.F/2)+1+View.Force.TimeOrder/2);
% h = bar3(View.Ft.Amplitude,0.6);
% zlim([0 100000]);
% for i = 1:numel(h)
%     zData = get(h(i),'ZData');
%     zData = repmat(max(zData,[],2),1,4);
%     set(h(i),'CData',zData);
%     set(h(i),'FaceColor','flat');
% end
% % imagesc(View.Ft.Amplitude);
% set(gca,'xtick',1:5:View.Force.TimeOrder+1);
% set(gca,'xticklabel',-View.Force.TimeOrder/2:5:View.Force.TimeOrder/2);
% title(s(4),'Tangential Force Density fft','Fontsize',14);
% % xlabel('Time Harmonic','Fontsize',14);
% set(gca,'ytick',1:5:View.Force.SpaceOrder+1);
% set(gca,'yticklabel',0:5:View.Force.SpaceOrder);
% % ylabel('Space Harmonic','Fontsize',14);
% zlabel('\itP\rm/(N/m^2)','Fontsize',14);
% set(gca,'linewidth',1);
% grid on
% box on
% view(3);

%% 磁密信号FFT分析
[Row.B,Column.B]=size(TimeDomain.Br);
% 径向磁密分解
Fourier.Br.ComplexData=fft2(TimeDomain.Br,N.Space,N.Time);
Fourier.Br.Amplitude=4*abs(Fourier.Br.ComplexData)/(N.Space*N.Time);% 幅值转化
Fourier.Br.Amplitude(1,:)=Fourier.Br.Amplitude(1,:)/2;% 直流分量幅值转化
Fourier.Br.Amplitude(:,1)=Fourier.Br.Amplitude(:,1)/2;% 直流分量幅值转化
Fourier.Br.Amplitude=fftshift(Fourier.Br.Amplitude);
Fourier.Br.ComplexData=fftshift(Fourier.Br.ComplexData);
Fourier.Br.Phase=angle(Fourier.Br.ComplexData)*180/pi;% 相位求解
% % 切向磁密分解
% Fourier.Bt.ComplexData=fft2(TimeDomain.Bt,N.Space,N.Time);
% Fourier.Bt.Amplitude=4*abs(Fourier.Bt.ComplexData)/(N.Space*N.Time);% 幅值转化
% Fourier.Bt.Amplitude(1,:)=Fourier.Bt.Amplitude(1,:)/2;% 直流分量幅值转化
% Fourier.Bt.Amplitude(:,1)=Fourier.Bt.Amplitude(:,1)/2;% 直流分量幅值转化
% Fourier.Bt.Amplitude=fftshift(Fourier.Bt.Amplitude);
% Fourier.Bt.ComplexData=fftshift(Fourier.Bt.ComplexData);
% Fourier.Bt.Phase=angle(Fourier.Bt.ComplexData)*180/pi;% 相位求解
% 
% Fourier.FluxDensity.Frequency=(0:1:N.Time/2)*Fs.Time/N.Time;% 时间频率转化
% Fourier.FluxDensity.TimeOrder=Fourier.FluxDensity.Frequency/Fourier.FluxDensity.Frequency(2);% 时间阶次求解
% Fourier.FluxDensity.SpaceOrder=(0:1:N.Space/2)*Fs.Space/N.Space;
% Fourier.SpaceOrder=Fourier.FluxDensity.SpaceOrder/Fourier.FluxDensity.SpaceOrder(2);% 空间阶次求解

%% 图形绘制
% 磁密分析结果显示
set(0,'defaultfigurecolor','w')
figure(2)
View.FluxDensity.SpaceOrder=48;
View.FluxDensity.TimeOrder=60;

% 径向磁密时域图
s(1)=subplot(2,2,1);
surface(TimeDomain.Time,TimeDomain.ElecAngle,TimeDomain.Br);
shading interp;
axis([0 3/325 0 360 -1 1]);
title(s(1),'Radial flux density','Fontsize',14);
zlabel('B(T)','Fontsize',14);
set(gca,'xtick',0:0.002:0.008);
set(gca,'xticklabel',0:2:8);
set(gca,'ytick',0:40:360);
set(gca,'yticklabel',0:40:360);
view(3);
box on
grid on
% 
% % 切向磁密时域图
% s(2)=subplot(2,2,2);
% surface(TimeDomain.Time,TimeDomain.ElecAngle,TimeDomain.Bt);
% shading interp;
% axis([0 3/325 0 360 -0.4 0.4]);
% title(s(2),'Tangential flux density','Fontsize',14);
% zlabel('B(T)','Fontsize',14);
% set(gca,'xtick',0:0.002:0.008);
% set(gca,'xticklabel',0:2:8);
% set(gca,'ytick',0:40:360);
% set(gca,'yticklabel',0:40:360);
% view(3);
% box on
% grid on
% 
% 径向磁密幅频图
s(3)=subplot(2,2,3);
View.Br.Amplitude=Fourier.Br.Amplitude(floor(Row.B/2)+1:floor(Row.B/2)+1+View.FluxDensity.SpaceOrder,floor(Column.B/2)+1-View.FluxDensity.TimeOrder/2:floor(Column.B/2)+1+View.FluxDensity.TimeOrder/2);
h = bar3(View.Br.Amplitude,0.6);
zlim([0 1.5]);
for i = 1:numel(h)
    zData = get(h(i),'ZData');
    zData = repmat(max(zData,[],2),1,4);
    set(h(i),'CData',zData);
    set(h(i),'FaceColor','flat');
end
% imagesc(View.Br.Amplitude);
set(gca,'xtick',1:5:View.FluxDensity.TimeOrder+1);
set(gca,'xticklabel',-View.FluxDensity.TimeOrder/2:5:View.FluxDensity.TimeOrder/2);
title(s(3),'Ridial Flux Density fft','Fontsize',14);
% xlabel('Time Harmonic','Fontsize',14);
set(gca,'ytick',1:5:View.FluxDensity.SpaceOrder+1);
set(gca,'yticklabel',0:5:View.FluxDensity.SpaceOrder);
% ylabel('Space Harmonic','Fontsize',14);
zlabel('B(T)','Fontsize',14);
set(gca,'linewidth',1);
grid on
box on
view(3);
% 
% % 切向磁密幅频图
% s(4)=subplot(2,2,4);
% View.Bt.Amplitude=Fourier.Bt.Amplitude(floor(Row.B/2)+1:floor(Row.B/2)+1+View.FluxDensity.SpaceOrder,floor(Column.B/2)+1-View.FluxDensity.TimeOrder/2:floor(Column.B/2)+1+View.FluxDensity.TimeOrder/2);
% h = bar3(View.Bt.Amplitude,0.6);
% zlim([0 0.5]);
% for i = 1:numel(h)
%     zData = get(h(i),'ZData');
%     zData = repmat(max(zData,[],2),1,4);
%     set(h(i),'CData',zData);
%     set(h(i),'FaceColor','flat');
% end
% % imagesc(View.Br.Amplitude);
% set(gca,'xtick',1:5:View.FluxDensity.TimeOrder+1);
% set(gca,'xticklabel',-View.FluxDensity.TimeOrder/2:5:View.FluxDensity.TimeOrder/2);
% title(s(3),'Tangential Flux Density fft','Fontsize',14);
% % xlabel('Time Harmonic','Fontsize',14);
% set(gca,'ytick',1:5:View.FluxDensity.SpaceOrder+1);
% set(gca,'yticklabel',0:5:View.FluxDensity.SpaceOrder);
% % ylabel('Space Harmonic','Fontsize',14);
% zlabel('B(T)','Fontsize',14);
% set(gca,'linewidth',1);
% grid on
% box on
% view(3);
end