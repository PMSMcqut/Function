%function [Fourier]=FFT_fun(data,time,space,f0,fmax,type,method)
clear;clc;
close all
load('OriginData.txt')
% load('FFT_signal.mat')
[r c]=size(OriginData);
[Fourier]=FFT_fun(OriginData(1:r,2:c),OriginData(1:r,1),1,50,18000,'1D','fun');
% ----- ͼ���������� -----------------------------------------------
fig=gcf;
NameArrayGcf={'Name','color','position','Units'};
VlaueArrayGcf={'Fourier Transfer','white',[500,100,800,500],'centimeters'};
set(fig,NameArrayGcf,VlaueArrayGcf);
% ------ Fourier Transfer -----------------------------------------------
subplot(1,2,1)
% ----- ͼ����������--------------------------------------------------
h1=plot(OriginData(1:r,1),OriginData(1:r,2:c));
NameArrayFig={'LineWidth','Marker','LineStyle','color'};
ValueArrayFig={1.5,'none','-','r'};% ����������ͼ���е��м�����ͼ�Ķ�������������ÿ�������м�������
set(h1,NameArrayFig,ValueArrayFig)
% ----- �������������� ----------------------------------------------------
ax=gca;
NameArrayAx={'FontSize','FontName','LineWidth','DataAspectRatio','xlim','xtick','xticklabel','ylim'};
ValueArrayAx={13.5,'Times New Roman',1,[1 1000 1],[0,0.06],0:0.01:0.06,0:0.01:0.06,[-20,20]};
set(ax,NameArrayAx,ValueArrayAx);
% ----- �������ǩ���� ----------------------------------------------------
tx=xlabel('');
NameArrayTx={'String','FontSize'};
ValueArrayTx={'\it time',14};
set(tx,NameArrayTx,ValueArrayTx);
ty=ylabel('');
NameArrayTy={'String','FontSize'};
ValueArrayTy={'\it Current',14};
set(ty,NameArrayTy,ValueArrayTy);
% title('1500 rpm �����߷����Ƽ�г������','FontName','YaHei Monaco Hybird','FontSize',11)
% ----- ͼ������ ---------------------------------------------------------
% lgd=legend;
% NameArrayLgd={'String','FontSize','Location'};
% ValueArrayLgd={{'Phase A','Phase D'},13.5,'northeast'};
% set(lgd,NameArrayLgd,ValueArrayLgd);


% ------ Fourier Transfer -----------------------------------------------
subplot(1,2,2)
% ----- ͼ����������--------------------------------------------------
h2=bar(Fourier.P.Order(1:400,:),Fourier.P.Amplitude(1:400,:));
NameArrayFig={'FaceColor','BarWidth'};
ValueArrayFig={'r',0.7};% ����������ͼ���е��м�����ͼ�Ķ�������������ÿ�������м�������
set(h2,NameArrayFig,ValueArrayFig)
% h=stem(x,y);
% NameArrayFig={'LineWidth','Marker','LineStyle','color'};
% ValueArrayFig={1.5,'*','none','k';
%     1.5,'x','-','r'};% ����������ͼ���е��м�����ͼ�Ķ�������������ÿ�������м�������
% set(h,NameArrayFig,ValueArrayFig)
% ----- �������������� ----------------------------------------------------
ax=gca;
NameArrayAx={'FontSize','FontName','LineWidth','DataAspectRatio','xlim','xtick','xticklabel','ylim'};
ValueArrayAx={13.5,'Times New Roman',1,[1 0.001 1],[-1,400],0:40:400,0:40:400,[0,0.4]};
set(ax,NameArrayAx,ValueArrayAx);
% ----- �������ǩ���� ----------------------------------------------------
tx=xlabel('');
NameArrayTx={'String','FontSize'};
ValueArrayTx={'\it Harmonic',14};
set(tx,NameArrayTx,ValueArrayTx);
ty=ylabel('');
NameArrayTy={'String','FontSize'};
ValueArrayTy={'\it Current',14};
set(ty,NameArrayTy,ValueArrayTy);
% ----- ͼ������ ---------------------------------------------------------
% lgd=legend;
% NameArrayLgd={'String','FontSize','Location'};
% ValueArrayLgd={{'Phase A','Phase D'},13.5,'northeast'};
% set(lgd,NameArrayLgd,ValueArrayLgd);
%% THD calculation
THD=sqrt(sum((Fourier.P.Amplitude(3:end)/Fourier.P.Amplitude(2)).^2));


% scrsz = get(groot, 'ScreenSize');
% set(groot, 'defaultFigurePosition', [10 100 scrsz(3)/3 scrsz(4)/3]);    % ����ͼ���С
% set(groot, 'defaultFigurePaperPositionMode', 'auto');                     % ����figure��ӡ��С�͵�ǰͼ���Сһ��
% set(groot, 'defaultFigureColor', 'white');                                 % ����figure��ɫ
% set(groot, 'defaultFigureNumberTitle', 'off');                             % ����figure���ⲻ��ʾĬ��figure����
% set(groot, 'defaultLineLineWidth', 1.5);                                    % ��������ͼ�����������Ϊ1
% set(groot, 'defaultAxesLineWidth', 1.5);                                    % ����ͼ��������ϸ
% set(groot, 'defaultAxesFontWeight', 'bold');                               % ����ͼ��������ϸ
% set(groot, 'defaultAxesFontName', 'times new roman');                     % ��������ͼ��������Ϊtimes new roman
% set(groot, 'defaultAxesFontSize', 14);                                     % �����figure������������Ĵ�С����ȥ����
% set(groot, 'defaultAxesXGrid', 'on');                                       % ����x��grid ��
% set(groot, 'defaultAxesYGrid', 'on');                                       % ����y��grid ��
% set(groot, 'defaultAxesBox', 'on');                                         % ����box on
% set(groot, 'defaultAxesGridLineStyle', ':');
% set(groot, 'defaultAxesGridAlpha', 0.35);