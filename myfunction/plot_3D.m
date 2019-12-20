% 3D-plot function: the plot stype includes surf, surfc, pcolor,contour,contourf,stem3 and bar3
function [fig,h,ax,label]=plot_3D(x,y,z,PlotStyle)
% ----- 图窗属性设置 -----------------------------------------------
fig=gcf;
NameArrayGcf={'Name','color','position','Units'};
VlaueArrayGcf={'3-D waveform of force density','white',[500,100,600,400],'centimeters'};
set(fig,NameArrayGcf,VlaueArrayGcf);
switch PlotStyle
    case 'surf'
        % ----- 图形属性设置-----------------------------------------
        h=surf(x,y,z);
        NameArrayGcf={'FaceAlpha','FaceColor','EdgeColor','LineStyle'};
        VlaueArrayGcf={1,'interp','interp','-'};
        set(h,NameArrayGcf,VlaueArrayGcf);
        % ----- 坐标轴属性设置 --------------------------------------
        ax=gca;
        NameArrayAx={'FontSize','FontName','LineWidth','Box','BoxStyle'};
        ValueArrayAx={13.5,'Times New Roman',1,'on','full'};
        set(ax,NameArrayAx,ValueArrayAx);
        % ----- 坐标轴标签设置 --------------------------------------
        label.tx=xlabel('');
        NameArrayTx={'String','FontSize','Position'};
        ValueArrayTx={'\it t \rm(s)',14,[0.01,-30,-100000]};
        set(label.tx,NameArrayTx,ValueArrayTx);
        label.ty=ylabel('');
        NameArrayTy={'String','FontSize','Position'};
        ValueArrayTy={'\it \theta \rm(\circ)',14,[-0.003,200,-80000]};
        set(label.ty,NameArrayTy,ValueArrayTy);
        label.tz=zlabel('');
        NameArrayTz={'String','FontSize'};
        ValueArrayTz={'\it F_r \rm(N/m^2)',14};
        set(label.tz,NameArrayTz,ValueArrayTz);
        label.cb=colorbar;
        NameArrayCb={'LineWidth','FontSize','Visible','Location'};
        ValueArrayCb={0.5,13,'off','north'};
        set(label.cb,NameArrayCb,ValueArrayCb);
        view(-45,20);
    case 'bar3'
        View.SpaceOrder=-12:2:12;
        View.TimeOrder=0:2:36;
        [angle_tf,angle_loc]=ismember(View.SpaceOrder,z.Real.SpaceOrder);
        [time_tf,time_loc]=ismember(View.TimeOrder,z.Real.TimeOrder);
        % ----- 图形属性设置-----------------------------------------
        h=bar3(z.Real.Amplitude(angle_loc,time_loc),0.5);
        zlim([0 max(max(z.Real.Amplitude))]);
        for i = 1:numel(h)
            zData = get(h(i),'ZData');
            zData = repmat(max(zData,[],2),1,4);
            set(h(i),'CData',zData);
            set(h(i),'FaceColor','flat');
        end
        % ----- 坐标轴属性设置 --------------------------------------
        ax=gca;
        NameArrayAx={'FontSize','FontName','LineWidth','Box','BoxStyle','xlim','xtick','xticklabel'...
            ,'ylim','ytick','yticklabel'};
        ValueArrayAx={13.5,'Times New Roman',1,'on','full',[0.5 length(time_loc)+0.5],1:2:length(time_loc),0:4:36,...
            [0.5 length(angle_loc)+0.5],1:2:length(angle_loc),-12:4:12};
        set(ax,NameArrayAx,ValueArrayAx);
        % ----- 坐标轴标签设置 --------------------------------------
        label.tx=xlabel('');
        NameArrayTx={'String','FontSize','Position'};
        ValueArrayTx={'\it f (\times f_0)',14,[12,12,-100000]};
        set(label.tx,NameArrayTx,ValueArrayTx);
        label.ty=ylabel('');
        NameArrayTy={'String','FontSize','Position'};
        ValueArrayTy={'\it r',14,[2,4,-80000]};
        set(label.ty,NameArrayTy,ValueArrayTy);
        label.tz=zlabel('');
        NameArrayTz={'String','FontSize'};
        ValueArrayTz={'\it F_r \rm(N/m^2)',14};
        set(label.tz,NameArrayTz,ValueArrayTz);
        label.cb=colorbar;
        NameArrayCb={'LineWidth','FontSize','Visible','Location'};
        ValueArrayCb={0.5,13,'off','north'};
        set(label.cb,NameArrayCb,ValueArrayCb);
        view(-45,20);
    case 'stem3'
        View.SpaceOrder=-12:2:12;
        View.TimeOrder=0:2:36;
        [angle_tf,angle_loc]=ismember(View.SpaceOrder,z.Real.SpaceOrder);
        [time_tf,time_loc]=ismember(View.TimeOrder,z.Real.TimeOrder);
        % ----- 图形属性设置-----------------------------------------
        h=stem3(z.Real.Amplitude(angle_loc,time_loc));
        NameArrayFig={'LineWidth','Marker','LineStyle','color'};
        ValueArrayFig={1.5,'x','-','r'};
        set(h,NameArrayFig,ValueArrayFig)
        % ----- 坐标轴属性设置 --------------------------------------
        ax=gca;
        NameArrayAx={'FontSize','FontName','LineWidth','Box','BoxStyle','xlim','xtick','xticklabel'...
            ,'ylim','ytick','yticklabel','YDir'};
        ValueArrayAx={13.5,'Times New Roman',1,'on','full',[0.5 length(time_loc)+0.5],1:2:length(time_loc),0:4:36,...
            [0.5 length(angle_loc)+0.5],1:2:length(angle_loc),-12:4:12,'reverse'};
        set(ax,NameArrayAx,ValueArrayAx);
        % ----- 坐标轴标签设置 --------------------------------------
        label.tx=xlabel('');
        NameArrayTx={'String','FontSize'};
        ValueArrayTx={'\it f (\times f_0)',14};
        set(label.tx,NameArrayTx,ValueArrayTx);
        label.ty=ylabel('');
        NameArrayTy={'String','FontSize'};
        ValueArrayTy={'\it r',14};
        set(label.ty,NameArrayTy,ValueArrayTy);
        label.tz=zlabel('');
        NameArrayTz={'String','FontSize'};
        ValueArrayTz={'\it F_r \rm(N/m^2)',14};
        set(label.tz,NameArrayTz,ValueArrayTz);
        label.cb=colorbar;
        NameArrayCb={'LineWidth','FontSize','Visible','Location'};
        ValueArrayCb={0.5,13,'off','north'};
        set(label.cb,NameArrayCb,ValueArrayCb);
        view(-45,20);
end
end