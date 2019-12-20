clear
clc
load('a1.mat');
load('a2.mat');
load('a4.mat');
save('alldata.mat');
load alldata;
step=[300 300 600];%转速间隔
%% a=1 n from 0 to 2700 rpm
n_final_a1=aa_nz_a1(:,[1:2700/step(1)+1]);
T_final_a1=aa_torque_a1(:,[1:2700/step(1)+1]);
Effi_final_a1=aa_efficiency_Pmax_a1(:,[1:2700/step(1)+1]);
%% a=2 n from 3000 to 10200 rpm
n_final_a2=aa_nz_a2(:,[3000/step(2)+1:10200/step(2)+1]);
T_final_a2=aa_torque_a2(:,[3000/step(2)+1:10200/step(2)+1]);
Effi_final_a2=aa_efficiency_Pmax_a2(:,[3000/step(2)+1:10200/step(2)+1]);
%% a=4 n from 10800 to 24000 rpm
n_final_a4=aa_nz_a4(:,[10800/step(3)+1:24000/step(3)+1]);
T_final_a4=aa_torque_a4(:,[10800/step(3)+1:24000/step(3)+1]);
Effi_final_a4=aa_efficiency_Pmax_a4(:,[10800/step(3)+1:24000/step(3)+1]);
%% 组合矩阵
n_final=[n_final_a1,n_final_a2,n_final_a4];
T_final=[T_final_a1,T_final_a2,T_final_a4];
Effi_final=[Effi_final_a1,Effi_final_a2,Effi_final_a4];
%% 绘制效率MAP
figure(1);
[C,h]=contour(n_final, T_final, Effi_final,50);
axis([0 24000 4 350]);%设置坐标轴刻度
colormap(jet);
hold on 
%% 绘制T-n曲线
[T_max T_index]=max(T_final,[],1);
[n_max n_index]=max(n_final,[],1);
plot(n_max,T_max,'b','Linewidth',2);
figure(2);
gca=pcolor(n_final, T_final, Effi_final);
axis([0 24000 4 350]);
colorbar;
colormap(jet);
hold on 
%% 绘制T-n曲线
[T_max T_index]=max(T_final,[],1);
[n_max n_index]=max(n_final,[],1);
plot(n_max,T_max,'b','Linewidth',2);
figure(3);
[C,h]=contourf(n_final, T_final, Effi_final,50);
axis([0 24000 5 350]);
colorbar;
colormap(jet);%设置图形颜色
clabel(C,h,'manual');%画等高线，条数用户自定义
hold on 
%% 绘制T-n曲线
[T_max T_index]=max(T_final,[],1);
[n_max n_index]=max(n_final,[],1);
plot(n_max,T_max,'b','Linewidth',2);
saveas(gcf,'final_efficiency_map');
%% 存储数据
save('final.mat');