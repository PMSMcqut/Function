
clear;
clc;
[NUM,TXT,RAW]=xlsread('iron loss-FEA-6000_3');
Id=NUM(:,1)*0.1;
Iq=NUM(:,2)*0.1;
xdata(1,:)=(Id.^2)';
xdata(2,:)=(Iq.^2)';
xdata(3,:)=Id;
xdata(4,:)=(Id.*Iq)';
xdata(5,:)=Iq;
xdata(6,:)=(Id.^3);
xdata(7,:)=(Iq.^3);
xdata(8,:)=(Id.^2.*Iq);
xdata(9,:)=(Iq.^2.*Id);

% xdata(4,1:49)=1;
% X=(xdata*xdata')^-1*xdata*NUM(:,[4,6]);

f=@(x,xda)x(1)*xda(1,:)+x(2)*xda(2,:)+x(3)*xda(3,:)+x(4)*xda(4,:)+x(5)*xda(5,:)+x(6)*xda(6,:)+x(7)*xda(7,:)+x(8)*xda(8,:)+x(9)*xda(9,:)+x(10);    %定义函数形式
%要拟合数据
ydata=NUM(:,3)';
a=[0,0,0,0,0,0,0,0,0,0];  %拟合初值
[x,resnorm]=lsqcurvefit(f,a,xdata,ydata)    %拟合
y=x(1)*xdata(1,:)+x(2)*xdata(2,:)+x(3)*xdata(3,:)+x(4)*xdata(4,:)+x(5)*xdata(5,:)+x(6)*xdata(6,:)+x(7)*xdata(7,:)+x(8)*xdata(8,:)+x(9)*xdata(9,:)+x(10);
bbb=y-ydata;
bbbb=abs(bbb)./ydata;
scatter(1:43,bbbb,3);
title('第6列');