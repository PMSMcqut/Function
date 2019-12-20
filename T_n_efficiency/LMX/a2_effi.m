%% Program Help Yang Lu 15,Oct.,2017  All rights reserved
% 1. 程序主要是将excel中转速，转矩，效率重新排列进行效率MAP的绘制
% 2. 在生成转矩矩阵和效率矩阵时，为了是图形连续，采用了将没有数据的地方补齐当前转矩最大值。
% 3. 不同并联之路时，转速间隔step需要改动，坐标轴的设置需要改动
clear;clc;
a=2; %并联支路数
step=300;%转速间隔
data=xlsread('efficiency.xlsx',['a=',num2str(a)]);%读取excel数据
l=length(find(data(:,1)==min(data(:,1))));%转速长度
n_max=max(data(:,1));%找出最大转速
zero=zeros(20,1);
n1=0:step:n_max;
a_nz=repmat(n1,l,1);%生成转速矩阵
T=[];effi=[];%定义转矩和效率矩阵
for n=0:step:n_max
    if n>0
    %% 求解生成转矩矩阵
    T_n=flip(data(find(data(:,1)==n),2));
    Tu=linspace(max(T_n),max(T_n),l-length(T_n));
    T(:,n/step+1)=[T_n;Tu'];
    Tmax(n/step+1)=max(T(:,n/step+1));
    %% 求解生成效率矩阵
    effi_n=flip(data(find(data(:,1)==n),3));
    effiu=linspace(max(effi_n),max(effi_n),l-length(effi_n));
    effi(:,n/step+1)=[effi_n;effiu'];
    else
    end
end
%% 绘图程序
aa_nz_a2=a_nz;aa_torque_a2=T;aa_efficiency_Pmax_a2=effi;%为了画图方便重置变量名称
%% 绘制T-n曲线
figure(1)
plot(n1,Tmax);
%% 绘制效率MAP
figure(1);
[C,h]=contour(aa_nz_a2, aa_torque_a2, aa_efficiency_Pmax_a2,50);
axis([0 12000 8 180]);%设置坐标轴刻度
colormap(jet);
hold on
plot(n1,Tmax,'b','Linewidth',2);
figure(2);
gca=pcolor(aa_nz_a2, aa_torque_a2, aa_efficiency_Pmax_a2);
axis([0 12000 8 180]);
colorbar;
colormap(jet);
hold on
plot(n1,Tmax,'b','Linewidth',2);
figure(3);
[C,h]=contourf(aa_nz_a2, aa_torque_a2, aa_efficiency_Pmax_a2,50);
axis([0 12000 8 180]);
colorbar;
colormap(jet);%设置图形颜色
clabel(C,h,'manual');%画等高线，条数用户自定义
hold on
plot(n1,Tmax,'b','Linewidth',2);
saveas(gcf,'efficiency_map_a2');
%% 存储数据
save('a2.mat','aa_nz_a2','aa_torque_a2','aa_efficiency_Pmax_a2');%存储数据，记得不同的并联之路要修改文件名
    
