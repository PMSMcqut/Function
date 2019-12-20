%function [Fourier]=FFT_fun(data,time,space,f0,fmax,type,method)
clear;clc;
close all
load('OriginData.txt')
% load('FFT_signal.mat')
[r c]=size(OriginData);
[Fourier]=FFT_fun(OriginData(1:r,2:c),OriginData(1:r,1),1,50,18000,'1D','fun');
% ----- 图窗属性设置 -----------------------------------------------
fig=gcf;
NameArrayGcf={'Name','color','position','Units'};
VlaueArrayGcf={'Fourier Transfer','white',[500,100,800,500],'centimeters'};
set(fig,NameArrayGcf,VlaueArrayGcf);
% ------ Fourier Transfer -----------------------------------------------
subplot(1,2,1)
% ----- 图形属性设置--------------------------------------------------
h1=plot(OriginData(1:r,1),OriginData(1:r,2:c));
NameArrayFig={'LineWidth','Marker','LineStyle','color'};
ValueArrayFig={1.5,'none','-','r'};% 行数代表了图形中的有几个绘图的对象，列数代表了每个对象有几个属性
set(h1,NameArrayFig,ValueArrayFig)
% ----- 坐标轴属性设置 ----------------------------------------------------
ax=gca;
NameArrayAx={'FontSize','FontName','LineWidth','DataAspectRatio','xlim','xtick','xticklabel','ylim'};
ValueArrayAx={13.5,'Times New Roman',1,[1 1000 1],[0,0.06],0:0.01:0.06,0:0.01:0.06,[-20,20]};
set(ax,NameArrayAx,ValueArrayAx);
% ----- 坐标轴标签设置 ----------------------------------------------------
tx=xlabel('');
NameArrayTx={'String','FontSize'};
ValueArrayTx={'\it time',14};
set(tx,NameArrayTx,ValueArrayTx);
ty=ylabel('');
NameArrayTy={'String','FontSize'};
ValueArrayTy={'\it Current',14};
set(ty,NameArrayTy,ValueArrayTy);
% title('1500 rpm 空载线反电势及谐波分析','FontName','YaHei Monaco Hybird','FontSize',11)
% ----- 图例设置 ---------------------------------------------------------
% lgd=legend;
% NameArrayLgd={'String','FontSize','Location'};
% ValueArrayLgd={{'Phase A','Phase D'},13.5,'northeast'};
% set(lgd,NameArrayLgd,ValueArrayLgd);


% ------ Fourier Transfer -----------------------------------------------
subplot(1,2,2)
% ----- 图形属性设置--------------------------------------------------
h2=bar(Fourier.P.Order(1:400,:),Fourier.P.Amplitude(1:400,:));
NameArrayFig={'FaceColor','BarWidth'};
ValueArrayFig={'r',0.7};% 行数代表了图形中的有几个绘图的对象，列数代表了每个对象有几个属性
set(h2,NameArrayFig,ValueArrayFig)
% h=stem(x,y);
% NameArrayFig={'LineWidth','Marker','LineStyle','color'};
% ValueArrayFig={1.5,'*','none','k';
%     1.5,'x','-','r'};% 行数代表了图形中的有几个绘图的对象，列数代表了每个对象有几个属性
% set(h,NameArrayFig,ValueArrayFig)
% ----- 坐标轴属性设置 ----------------------------------------------------
ax=gca;
NameArrayAx={'FontSize','FontName','LineWidth','DataAspectRatio','xlim','xtick','xticklabel','ylim'};
ValueArrayAx={13.5,'Times New Roman',1,[1 0.001 1],[-1,400],0:40:400,0:40:400,[0,0.4]};
set(ax,NameArrayAx,ValueArrayAx);
% ----- 坐标轴标签设置 ----------------------------------------------------
tx=xlabel('');
NameArrayTx={'String','FontSize'};
ValueArrayTx={'\it Harmonic',14};
set(tx,NameArrayTx,ValueArrayTx);
ty=ylabel('');
NameArrayTy={'String','FontSize'};
ValueArrayTy={'\it Current',14};
set(ty,NameArrayTy,ValueArrayTy);
% ----- 图例设置 ---------------------------------------------------------
% lgd=legend;
% NameArrayLgd={'String','FontSize','Location'};
% ValueArrayLgd={{'Phase A','Phase D'},13.5,'northeast'};
% set(lgd,NameArrayLgd,ValueArrayLgd);
%% THD calculation
THD=sqrt(sum((Fourier.P.Amplitude(3:end)/Fourier.P.Amplitude(2)).^2));


% scrsz = get(groot, 'ScreenSize');
% set(groot, 'defaultFigurePosition', [10 100 scrsz(3)/3 scrsz(4)/3]);    % 设置图框大小
% set(groot, 'defaultFigurePaperPositionMode', 'auto');                     % 设置figure打印大小和当前图框大小一致
% set(groot, 'defaultFigureColor', 'white');                                 % 设置figure底色
% set(groot, 'defaultFigureNumberTitle', 'off');                             % 设置figure标题不显示默认figure数字
% set(groot, 'defaultLineLineWidth', 1.5);                                    % 设置所绘图形中线条宽度为1
% set(groot, 'defaultAxesLineWidth', 1.5);                                    % 设置图框线条粗细
% set(groot, 'defaultAxesFontWeight', 'bold');                               % 设置图框线条粗细
% set(groot, 'defaultAxesFontName', 'times new roman');                     % 设置所绘图形中字体为times new roman
% set(groot, 'defaultAxesFontSize', 14);                                     % 会控制figure区内所有字体的大小，舍去不用
% set(groot, 'defaultAxesXGrid', 'on');                                       % 设置x轴grid 开
% set(groot, 'defaultAxesYGrid', 'on');                                       % 设置y轴grid 开
% set(groot, 'defaultAxesBox', 'on');                                         % 设置box on
% set(groot, 'defaultAxesGridLineStyle', ':');
% set(groot, 'defaultAxesGridAlpha', 0.35);