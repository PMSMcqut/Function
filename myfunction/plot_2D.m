% 2D-plot function: the plot stype includes line, stem and bar
function [fig,h,ax,tx,ty,lgd]=plot_2D(x,y,PlotStyle)
fig=gcf;
NameArrayGcf={'Name','color','position','Units'};
VlaueArrayGcf={'Spectrum of Current','white',[500,100,800,500],'centimeters'};
set(fig,NameArrayGcf,VlaueArrayGcf);
% ----- 图形属性设置--------------------------------------------------
switch PlotStyle
    case 'line'
        h=plot(x,y);
        NameArrayFig={'LineWidth','Marker','LineStyle','color'};
        ValueArrayFig={1.5,'none','-.','k';
            1.5,'none','-','r'};% 行数代表了图形中的有几个绘图的对象，列数代表了每个对象有几个属性
        set(h,NameArrayFig,ValueArrayFig)
    case 'stem'
        h=stem(x,y);
        NameArrayFig={'LineWidth','Marker','LineStyle','color'};
        ValueArrayFig={1.5,'*','none','k';
            1.5,'x','-','r'};% 行数代表了图形中的有几个绘图的对象，列数代表了每个对象有几个属性
        set(h,NameArrayFig,ValueArrayFig)
    case 'bar'
        h=bar(x,y);
        NameArrayFig={'FaceColor','BarWidth'};
        ValueArrayFig={'k',0.7;'r',0.7};% 行数代表了图形中的有几个绘图的对象，列数代表了每个对象有几个属性
        set(h,NameArrayFig,ValueArrayFig)
end
% ----- 坐标轴属性设置 ----------------------------------------------------
ax=gca;
NameArrayAx={'FontSize','FontName','LineWidth','DataAspectRatio','xlim','xtick','xticklabel','ylim'};
ValueArrayAx={13.5,'Times New Roman',1,[3000 1 1],[0,30],0:3:30,0:30,[0,0.004]};
set(ax,NameArrayAx,ValueArrayAx);
% ----- 坐标轴标签设置 ----------------------------------------------------
tx=xlabel('');
NameArrayTx={'String','FontSize'};
ValueArrayTx={'\it f \rm(Hz)',14};
set(tx,NameArrayTx,ValueArrayTx);
ty=ylabel('');
NameArrayTy={'String','FontSize'};
ValueArrayTy={'\it I \rm(A)',14};
set(ty,NameArrayTy,ValueArrayTy);
% ----- 图例设置 ---------------------------------------------------------
lgd=legend;
NameArrayLgd={'String','FontSize','Location'};
ValueArrayLgd={{'3phase','6phase'},13.5,'northwest'};
set(lgd,NameArrayLgd,ValueArrayLgd);
end