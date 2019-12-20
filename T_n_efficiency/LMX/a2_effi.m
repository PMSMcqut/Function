%% Program Help Yang Lu 15,Oct.,2017  All rights reserved
% 1. ������Ҫ�ǽ�excel��ת�٣�ת�أ�Ч���������н���Ч��MAP�Ļ���
% 2. ������ת�ؾ����Ч�ʾ���ʱ��Ϊ����ͼ�������������˽�û�����ݵĵط����뵱ǰת�����ֵ��
% 3. ��ͬ����֮·ʱ��ת�ټ��step��Ҫ�Ķ����������������Ҫ�Ķ�
clear;clc;
a=2; %����֧·��
step=300;%ת�ټ��
data=xlsread('efficiency.xlsx',['a=',num2str(a)]);%��ȡexcel����
l=length(find(data(:,1)==min(data(:,1))));%ת�ٳ���
n_max=max(data(:,1));%�ҳ����ת��
zero=zeros(20,1);
n1=0:step:n_max;
a_nz=repmat(n1,l,1);%����ת�پ���
T=[];effi=[];%����ת�غ�Ч�ʾ���
for n=0:step:n_max
    if n>0
    %% �������ת�ؾ���
    T_n=flip(data(find(data(:,1)==n),2));
    Tu=linspace(max(T_n),max(T_n),l-length(T_n));
    T(:,n/step+1)=[T_n;Tu'];
    Tmax(n/step+1)=max(T(:,n/step+1));
    %% �������Ч�ʾ���
    effi_n=flip(data(find(data(:,1)==n),3));
    effiu=linspace(max(effi_n),max(effi_n),l-length(effi_n));
    effi(:,n/step+1)=[effi_n;effiu'];
    else
    end
end
%% ��ͼ����
aa_nz_a2=a_nz;aa_torque_a2=T;aa_efficiency_Pmax_a2=effi;%Ϊ�˻�ͼ�������ñ�������
%% ����T-n����
figure(1)
plot(n1,Tmax);
%% ����Ч��MAP
figure(1);
[C,h]=contour(aa_nz_a2, aa_torque_a2, aa_efficiency_Pmax_a2,50);
axis([0 12000 8 180]);%����������̶�
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
colormap(jet);%����ͼ����ɫ
clabel(C,h,'manual');%���ȸ��ߣ������û��Զ���
hold on
plot(n1,Tmax,'b','Linewidth',2);
saveas(gcf,'efficiency_map_a2');
%% �洢����
save('a2.mat','aa_nz_a2','aa_torque_a2','aa_efficiency_Pmax_a2');%�洢���ݣ��ǵò�ͬ�Ĳ���֮·Ҫ�޸��ļ���
    
