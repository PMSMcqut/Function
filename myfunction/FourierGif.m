clear
h=figure(1);
filename='test.gif';
fps=30;
time=1:1*2*3*4*3;
cmap1=lines(8);
radius1=0.3./[-1 2 -3 4];
theta1=[1;2;3;4]*time/time(end)*2*pi;
locX1=radius1'*ones(1,numel(time)).*cos(theta1);
locY1=radius1'*ones(1,numel(time)).*sin(theta1);
locX2=zeros(5,numel(time));
locY2=zeros(5,numel(time));
for i=1:4
    locX2(i+1,:)=locX1(i,:)+locX2(i,:);
    locY2(i+1,:)=locY1(i,:)+locY2(i,:);
end
str1={'$$\frac{2sin\theta}{-\pi}$$','$$\frac{2sin2\theta}{2\pi}$$','$$\frac{2sin3\theta}{-3\pi}$$','$$\frac{2sin4\theta}{4\pi}$$'};
for i=time
    clf
    axes('position',[0 0 1 1],'ycolor','white','xcolor','white');
    hold on
    for j=1:4
        text(-1.2,-j,str1{j},'interpreter','latex','fontsize',16)
        plot(radius1(j)*cos((0:100)/100*2*pi),radius1(j)*sin((0:100)/100*2*pi)-j,'color',cmap1(j,:),'linewidth',2)
        plot([0 locX1(j,i)],[0 locY1(j,i)]-j,'color',cmap1(j,:),'linewidth',2)
        for k=1:j
            plot(radius1(k)*cos((0:100)/100*2*pi)+locX2(k,i)+1,radius1(k)*sin((0:100)/100*2*pi)+locY2(k,i)-j,'color',cmap1(k,:),'linewidth',2)
            plot(locX2(k+[0 1],i)+1,locY2(k+[0 1],i)-j,'color',cmap1(k,:),'linewidth',2)
        end
        plot(2+(time-1)/time(end),locY2(j+1,[i:-1:1 time(end):-1:i+1])-j,'color',cmap1(j,:),'linewidth',2)
        plot([locX2(j+1,i)+1 2],locY2(j+1,i)-j+[0 0],'color',cmap1(j,:),'linewidth',1)
        plot(2+[-0.1 0 -0.1],locY2(j+1,i)-j+[0.05 0 -0.05],'color',cmap1(j,:),'linewidth',1)
    end
    axis([-1.5 3.5 -4.6 -0.6])
    drawnow
    % Capture the plot as an image
    frame=getframe(h);
    im=frame2im(frame);
    [imind,cm]=rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1/fps);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/fps);
    end
end