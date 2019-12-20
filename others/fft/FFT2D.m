
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
%% ���ݶ�ȡ�봦��
p=3;
Qs=36;
Nt=gcd(Qs,p);
n=1500;%speed
f=p*n/60;%fundamental electrical frequency
T=1/f;%period
% OriginData.Br=xlsread('E:\Project\AFPM-Vibration\Force\DDW_bz_load.xlsx');
% OriginData.Bt=xlsread('E:\Project\AFPM-Vibration\Force\DDW_bt_load.xlsx');
OriginData.Br=load('E:\Project\AFPM-Vibration\Force\DDW_bz_load.xlsx');
OriginData.Bt=load('E:\Project\AFPM-Vibration\Force\DDW_bt_load.xlsx');
[Row.B,Column.B]=size(OriginData.Br);
TimeDomain.Br=OriginData.Br(:,2:Column.B);
TimeDomain.Bt=OriginData.Bt(:,2:Column.B);

TimeDomain.ElecAngle=OriginData.Br(:,1)*Nt;% �ռ��������
TimeDomain.Time=0:T/(Column.B-2):T;% ʱ���������
Fs.Time=1/(1/(p*n/60)/(Column.B-1));% ʱ�����Ƶ��
N.Time=length(TimeDomain.Time);% ʱ���������
Fs.Space=1/(TimeDomain.ElecAngle(2)-TimeDomain.ElecAngle(1));% �ռ����Ƶ��
N.Space=length(TimeDomain.ElecAngle);
[TimeDomain,Row,Column,Fourier]=ForceFluxCalculation(TimeDomain,N,Fs);
Fr_brbt_load=Fourier.Fr;
save('Fr_brbt_load.mat','Fr_brbt_load');