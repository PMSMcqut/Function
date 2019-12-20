clear all
clc
load in_power.dat;
%in_power中的数据必须是6列，而且从左到右依次为id,iq,Lamda_dm,Lamda_qm,Lamda_ds,Lamda_qs,
%而且Lamda_d可能是负值，并且六列中id是从小到大递增的
[N,N_row]=size(in_power);
%
id_ini=0;
iq_ini=0;
iq_max=120;
i_diff=20;
%以上参数每次需修改


d=20;%很神奇的一点，至今不知为什么，d不能比i_diff大
z=1;
eend=0;
data_num=N+1;
%id由小到大递增，将同一个id值下的不同iq继续细化离散
for i=1:(N-1)
   if in_power(i,1)==in_power(i+1,1)
       for m=0:(d-1)
           ii=d*(i-z)+m+z;
           ele(ii,1)=in_power(i,1);
           ele(ii,2)=in_power(i,2)+m*(in_power(i+1,2)-in_power(i,2))/d;
           ele(ii,3)=in_power(i,3)+m*(in_power(i+1,3)-in_power(i,3))/d;
           ele(ii,4)=in_power(i,4)+m*(in_power(i+1,4)-in_power(i,4))/d;
           ele(ii,5)=in_power(i,5)+m*(in_power(i+1,5)-in_power(i,5))/d;
           ele(ii,6)=in_power(i,6)+m*(in_power(i+1,6)-in_power(i,6))/d;
           
       end
       eend=eend+1;
   else
       ele((d*eend+z),1)=in_power(i,1);
       ele((d*eend+z),2)=in_power(i,2);
       ele((d*eend+z),3)=in_power(i,3);
       ele((d*eend+z),4)=in_power(i,4);
       ele((d*eend+z),5)=in_power(i,5);
       ele((d*eend+z),6)=in_power(i,6);
       z=z+1;
   end

end

ele(ii+1,1)=in_power(i+1,1);
ele(ii+1,2)=in_power(i+1,2);
ele(ii+1,3)=in_power(i+1,3);
ele(ii+1,4)=in_power(i+1,4);
ele(ii+1,5)=in_power(i+1,5);
ele(ii+1,6)=in_power(i+1,6);

%将上一步离散好的数据按照iq由小到大的顺序重新排列
[line_ele,row_ele]=size(ele);
j=1;
iq=iq_ini;
while iq<(iq_max+i_diff/d)
    for i=1:line_ele
        %if ele(i,1)==0
         %   break
        %end
        if ele(i,2)==iq
            ele1(j,1)=ele(i,1);
            ele1(j,2)=ele(i,2);
            ele1(j,3)=ele(i,3);
            ele1(j,4)=ele(i,4);
            ele1(j,5)=ele(i,5);
            ele1(j,6)=ele(i,6);
            j=j+1;
        end
    
    end
    iq=iq+i_diff/d;
end    

z=1;
eend=0;
%将重新排列的数据，在iq相同的情况下，插值得到进一步扩展的数据
for i=1:(line_ele-1)
   %if ele1(i,1)==0
    %   break
   %end
   if ele1(i,2)==ele1(i+1,2)
       for m=0:(d-1)
           ii=d*(i-z)+m+z;
           i_d(ii)=ele1(i,1)+m*(ele1(i+1,1)-ele1(i,1))/d;
           i_q(ii)=ele1(i,2);
           Lamda_dm(ii)=ele1(i,3)+m*(ele1(i+1,3)-ele1(i,3))/d;
           Lamda_qm(ii)=ele1(i,4)+m*(ele1(i+1,4)-ele1(i,4))/d;
           Lamda_ds(ii)=ele1(i,5)+m*(ele1(i+1,5)-ele1(i,5))/d;
           Lamda_qs(ii)=ele1(i,6)+m*(ele1(i+1,6)-ele1(i,6))/d; 
           
            ele2(ii,1)=i_d(ii);
            ele2(ii,2)=i_q(ii);
            ele2(ii,3)=Lamda_dm(ii);
            ele2(ii,4)=Lamda_qm(ii);
            ele2(ii,5)=Lamda_ds(ii);
            ele2(ii,6)=Lamda_qs(ii);
           
       end
       eend=eend+1;
   else
       i_d(d*eend+z)=ele1(i,1);
       i_q(d*eend+z)=ele1(i,2);
       Lamda_dm(d*eend+z)=ele1(i,3);
       Lamda_qm(d*eend+z)=ele1(i,4);
       Lamda_ds(d*eend+z)=ele1(i,5);
       Lamda_qs(d*eend+z)=ele1(i,6);
       
            ele2(d*eend+z,1)=i_d(d*eend+z);
            ele2(d*eend+z,2)=i_q(d*eend+z);
            ele2(d*eend+z,3)=Lamda_dm(d*eend+z);
            ele2(d*eend+z,4)=Lamda_qm(d*eend+z);
            ele2(d*eend+z,5)=Lamda_ds(d*eend+z);
            ele2(d*eend+z,6)=Lamda_qs(d*eend+z);
       z=z+1;
   end

end
       i_d(ii+1)=ele1(i+1,1);
       i_q(ii+1)=ele1(i+1,2);
       Lamda_dm(ii+1)=ele1(i+1,3);
       Lamda_qm(ii+1)=ele1(i+1,4);
       Lamda_ds(ii+1)=ele1(i+1,5);
       Lamda_qs(ii+1)=ele1(i+1,6);   
            ele2(ii+1,1)=i_d(ii+1);
            ele2(ii+1,2)=i_q(ii+1);
            ele2(ii+1,3)=Lamda_dm(ii+1);
            ele2(ii+1,4)=Lamda_qm(ii+1);
            ele2(ii+1,5)=Lamda_ds(ii+1);
            ele2(ii+1,6)=Lamda_qs(ii+1);   
 
[line_ele2,row_ele2]=size(ele2); 

%以下为电机参数输入及特性曲线计算绘图
mode=1;%mode=0:motoring;else:generatoring %%%%mode=0是发电机模型
if mode==0
n0=0;%n0-初始转速
n_max=12000;
n_step=100;%转速间隔100转
nn=n_max/n_step+1-n0/n_step;

p=2;%极对数
i_max=11;
v_max=310*0.85/1.732;%相电压？
Rs=0.7561;%
 
 
for kh=1:nn
    n=n0+n_step*(kh-1);
    nz(kh)=n;
    f=n*p/60;%电角频率
    fz(kh)=f;
   
    w=2*pi*f;%电角速度
    wz(kh)=w;
    
    power_max=0/1000;%最大输出功率,起始值应设为0，后面的程序有一直在找最大功率
    ibb=10;%在in_power的数据里，id，iq被放大了ibb倍，所以在下面的计算里要把这个倍数除回来
for i=1:line_ele2
    id=ele2(i,1)/ibb;
    iq=ele2(i,2)/ibb;
   
    Ld_m=ele2(i,3);% m-应该指空载,L是指Lamda
   
    Lq_m=ele2(i,4);% s-应该指负载
   
    Ld_s=ele2(i,5);
    
    Lq_s=ele2(i,6);
  
    
       % if id==0
       %     break
       % else
            Iref=sqrt(id^2+iq^2);
            
            if Iref<i_max
                 Vdm=(-1)*w*Lq_m; %正比于角速度的运动电势
                 Vds=(-1)*w*Lq_s;  %负载运动电势
                 Vqm=(1)*w*(1)*Ld_m;  
                 Vqs=(1)*w*(1)*Ld_s;
                 Vdmz(i,kh)=Vdm;
                 Vdsz(i,kh)=Vds;  
                 Vqmz(i,kh)=Vqm;  
                 Vqsz(i,kh)=Vqs;
                
                 Vdm_t=(+1)*w*Lq_m-Rs*id; %%%%改了最后一个负号
                 Vds_t=(+1)*w*Lq_s-Rs*id;  %%%%改了最后一个负号
                 Vqm_t=(-1)*w*(1)*Ld_m-Rs*iq;  
                 Vqs_t=(-1)*w*(1)*Ld_s-Rs*iq;
                 Vdm_tz(i,kh)=Vdm_t;
                 Vds_tz(i,kh)=Vds_t;
                 Vqm_tz(i,kh)=Vqm_t;
                 Vqs_tz(i,kh)=Vqs_t;
                 
                 power_m=-3/2*((1)*Vdm*id+Vqm*iq)/1000;%单位是kW  %%%加了负号
                 power_s=-3/2*((1)*Vds*id+Vqs*iq)/1000;
                 power_im=3/2*((-1)*Vqm*id-Vdm*iq)/1000;
                 power_is=3/2*((-1)*Vqs*id-Vds*iq)/1000;
                 power_factor_m= power_m/sqrt(power_m^2+power_im^2);
                 power_factor_s= power_s/sqrt(power_s^2+power_is^2);
                 power_mz(i,kh)=power_m/1000;
                 power_sz(i,kh)=power_s/1000;
                 power_imz(i,kh)=power_im/1000;
                 power_isz(i,kh)=power_is/1000;
                 power_factor_mz(i,kh)=power_factor_m; 
                 power_factor_sz(i,kh)=power_factor_s; 
                 
                 efficiency_m=power_m/(power_m+1.5*(id^2+iq^2)*Rs/1000 + 2.000*(n/n_max)^2)*100;%2000是什么？
                 efficiency_s=power_s/(power_s+1.5*(id^2+iq^2)*Rs/1000 + 2.000*(n/n_max)^2)*100;
                 Vabs_m=sqrt(Vdm_t^2+Vqm_t^2);
                 Vabs_s=sqrt(Vds_t^2+Vqs_t^2);
                 efficiency_mz(i,kh)=efficiency_m;
                 efficiency_sz(i,kh)=efficiency_s;
                 Vabs_mz(i,kh)=Vabs_m;
                 Vabs_sz(i,kh)=Vabs_s;
                 
                 if Vabs_s<v_max
                     %if efficiency_s>45
                        %a_efficiency_s(i,kh)=efficiency_s;
                        %a_torque_s(i,kh)=power_s/(w/p)*1000;
                     %end
                     
                     if power_s>power_max
                        Id=id;
                        a_Id(kh)=Id;
                        Iq=iq;
                        a_Iq(kh)=Iq;
                        Imax=sqrt(Id^2+Iq^2);
                        a_Imax(kh)=Imax;
                        Vpeak=Vabs_s;
                        a_Vpeak(kh)=Vpeak;
                        Vq=Vqs;
                        a_Vq(kh)=Vq;
                        Vd=Vds;
                        a_Vd(kh)=Vd;
                        Vq_t=Vqs_t;
                        a_Vq_t(kh)=Vq_t;
                        Vd_t=Vds_t;
                        a_Vd_t(kh)=Vd_t;
                        a_Ld_s(kh)=Ld_s;
                        a_Lq_s(kh)=Lq_s;
                        power_max=power_s;
                        a_power_max(kh)=power_max;
                        a_power_max_im(kh)=power_is;
                        a_power_max_t(kh)=power_s+1.5*(id.^2+iq.^2)*Rs;
                        a_efficiency_Pmax(kh)=efficiency_s;
					    a_p_factor(kh)=power_factor_s;
				        a_torque(kh)=power_max/(w/p)*1000;
                        angle=acosd(-Id/Imax);
                        a_angle(kh)=angle;
                     end
                 end
                 
            end
        
        end
    end


else
n0=0;%n0-初始转速
n_max=3995;
n_step=100;%转速间隔100转
nn=n_max/n_step+1-n0/n_step;

p=2;%极对数
Rs=0.7561;%
i_max=11;
v_max=310*0.85/1.732;%相电压？
%v_max=117;
for kh=1:nn
    n=n0+n_step*(kh-1);
    nz(kh)=n;
    f=n*p/60;%电角频率
    fz(kh)=f;
    
    w=2*pi*f;%电角速度
    wz(kh)=w;
   
    power_max=0/1000;%最大输出功率,起始值应设为0，后面的程序有一直在找最大功率
for i=1:line_ele2
    id=ele2(i,1)/10;
    iq=ele2(i,2)/10;
   
    Ld_m=ele2(i,3);% m-应该指空载,L是指Lamda
   
    Lq_m=ele2(i,4);% s-应该指负载
   
    Ld_s=ele2(i,5);
    
    Lq_s=ele2(i,6);
  
    
       % if id==0
       %     break
       % else
            Iref=sqrt(id^2+iq^2);
            
            if Iref<i_max
                 Vdm=(-1)*w*Lq_m; %正比于角速度的运动电势
                 Vds=(-1)*w*Lq_s;  
                 Vqm=(1)*w*(1)*Ld_m;  
                 Vqs=(1)*w*(1)*Ld_s;
                 Vdmz(i,kh)=Vdm;
                 Vdsz(i,kh)=Vds;  
                 Vqmz(i,kh)=Vqm;  
                 Vqsz(i,kh)=Vqs;
                
                 Vdm_t=(-1)*w*Lq_m+Rs*id; 
                 Vds_t=(-1)*w*Lq_s+Rs*id;  
                 Vqm_t=(1)*w*(1)*Ld_m+Rs*iq;  
                 Vqs_t=(1)*w*(1)*Ld_s+Rs*iq;
                 Vdm_tz(i,kh)=Vdm_t;
                 Vds_tz(i,kh)=Vds_t;
                 Vqm_tz(i,kh)=Vqm_t;
                 Vqs_tz(i,kh)=Vqs_t;
                 
                 power_m=3/2*((1)*Vdm*id+Vqm*iq)/1000;%单位是kW   
                 power_s=3/2*((1)*Vds*id+Vqs*iq)/1000;   
                 power_im=3/2*((1)*Vqm*id+Vdm*iq)/1000;   
                 power_is=3/2*((1)*Vqs*id+Vds*iq)/1000;   
                 power_factor_m= power_m/sqrt(power_m^2+power_im^2);
                 power_factor_s= power_s/sqrt(power_s^2+power_is^2);
                 power_mz(i,kh)=power_m/1000;
                 power_sz(i,kh)=power_s/1000;
                 power_imz(i,kh)=power_im/1000;
                 power_isz(i,kh)=power_is/1000;
                 power_factor_mz(i,kh)=power_factor_m; 
                 power_factor_sz(i,kh)=power_factor_s; 
                 Vabs_m=sqrt(Vdm_t^2+Vqm_t^2);
                 Vabs_s=sqrt(Vds_t^2+Vqs_t^2);
                 Vabs_mz(i,kh)=Vabs_m;
                 Vabs_sz(i,kh)=Vabs_s;
                 n_ref=2000;
                 efficiency_m=(-power_m)/(-power_m+1.5*(id^2+iq^2)*Rs/1000 + 10.24/1000*(n/n_ref)^2)*100;%2000是什么？
                 
                 efficiency_s=(-power_s)/(-power_s+1.5*(id^2+iq^2)*Rs/1000 + 10.24/1000*(n/n_ref)^2)*100;
                 
                 
                 efficiency_mz(i,kh)=efficiency_m;
                 efficiency_sz(i,kh)=efficiency_s;
                 
                 
                 if Vabs_s<v_max
                     %if efficiency_s>45
                        %a_efficiency_s(i,kh)=efficiency_s;
                        %a_torque_s(i,kh)=power_s/(w/p)*1000;
                     %end
                     
                     if power_s<power_max
                        Id=id;
                        a_Id(kh)=Id;
                        Iq=iq;
                        a_Iq(kh)=Iq;
                        Imax=sqrt(Id^2+Iq^2);
                        a_Imax(kh)=Imax;
                        Vpeak=Vabs_s;
                        a_Vpeak(kh)=Vpeak;
                        Vq=Vqs;
                        a_Vq(kh)=Vq;
                        Vd=Vds;
                        a_Vd(kh)=Vd;
                        Vq_t=Vqs_t;
                        a_Vq_t(kh)=Vq_t;
                        Vd_t=Vds_t;
                        a_Vd_t(kh)=Vd_t;
                        a_Ld_s(kh)=Ld_s;
                        a_Lq_s(kh)=Lq_s;
                        power_max=power_s;
                        a_power_max(kh)=-power_max;  %
                        a_power_max_im(kh)=power_is;
                        a_power_max_t(kh)=-power_s-1.5*(id.^2+iq.^2)*Rs/1000;  %
                        a_efficiency_Pmax(kh)=efficiency_s;
					    a_p_factor(kh)=power_factor_s;
				        a_torque(kh)=-power_max/(w/p)*1000;  %
                        angle=acosd(-Id/Imax);
                        a_angle(kh)=angle;
                     end
                 end
                 
            end
        
        end
    end
 end


F1=figure(1);
ylabels{1} = 'Torque (N*m)';
ylabels{2} = 'Power (kW)';
ylabels{3} = 'Vmax (V)';
ylabels{4} = 'Angle (deg.）';
%Plot the first two lines with plotyy
[ax,hlines(1),hlines(2)] = plotyy(nz,a_torque,nz,a_power_max);
xlabel('n (rpm)','fontsize',16);
set(hlines(1),'linestyle','-','marker','o','LineWidth',2)
set(hlines(2),'linestyle','-','marker','p','LineWidth',2)
set(ax(1),'Fontsize',16,'ytick',[0:1:35]);
set(ax(2),'Fontsize',16,'ytick',[0:0.4:10]);
cfig = get(gcf,'color');
pos = [0.125 0.1 0.65 0.8];
offset = pos(3)/5.5;
%Reduce width of the two axes generated by plotyy
pos(1) = pos(1) + offset;
pos(3) = pos(3) - offset;
set(ax,'position',pos);
%Determine the position of the third/fourth axes
pos3 = [pos(1) pos(2) pos(3)+offset pos(4)];
pos4 = [pos(1) - offset pos(2) pos(3)+offset pos(4)];
%Determine the proper x-limits for the third and fourth axes
scale3 = pos3(3)/pos(3);
scale4 = pos4(3)/pos(3);
limx1 = get(ax(1),'xlim');
limx3 = [limx1(1) limx1(1)+scale3*(limx1(2)-limx1(1))];
limx4 = [limx1(2)-scale4*(limx1(2)-limx1(1)) limx1(2)];

%Create ax(3) & ax(4)
ax(3) = axes('Position',pos3,'box','off',...
   'Color','none','XColor',cfig,'YColor','r',...
   'xtick',[],'xlim',limx3,'yaxislocation','right');
ax(4) = axes('Position',pos4,'box','off',...
   'Color','none','XColor',cfig,'YColor','k',...
   'xtick',[],'xlim',limx4,'yaxislocation','left');

%Plot x3,y3,x4,y4
hlines(3) = line(nz,a_Vpeak,'Parent',ax(3));

hlines(4) = line(nz,a_angle,'Parent',ax(4));

set(hlines(3),'linestyle','-','marker','<','color','r','LineWidth',2)
set(hlines(4),'linestyle','-','marker','s','color','k','LineWidth',2)
set(ax(3),'Fontsize',16)
set(ax(4),'Fontsize',16)
%Set y-labels;
set(get(ax(1),'ylabel'),'string',ylabels{1});
set(get(ax(2),'ylabel'),'string',ylabels{2});
set(get(ax(3),'ylabel'),'string',ylabels{3});
set(get(ax(4),'ylabel'),'string',ylabels{4});
legend([hlines(1),hlines(2),hlines(3),hlines(4)],'Torque','Power','Vmax','Angle',2);
gca=legend([hlines(1),hlines(2),hlines(3),hlines(4)],'Torque','Power','Vmax','Angle',2);
set(gca,'position',[0.409722222982799 0.260989962698224 0.0546874992394198 0.0314987877344499]);
set( F1, 'Color','w');

%画效率—转矩-转速图

%删除转矩、效率全零行
%a_efficiency_s(all(a_efficiency_s==0,2),:)=[];
%a_torque_s(all(a_torque_s==0,2),:)=[];

%[line_efficiency,row_efficiency]=size(a_efficiency_s);
%a_nz=repmat(nz,line_efficiency,1);
%figure(2);
%[C,h]=contour(a_nz, a_torque_s, a_efficiency_s,100);


%figure(3);
%gca=pcolor(a_nz, a_torque_s, a_efficiency_s);
%set(gca, 'LineStyle','none');
%colorbar;

