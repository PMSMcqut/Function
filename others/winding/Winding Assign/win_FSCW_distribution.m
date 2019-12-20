% calculation of three-phase winding distribution and mmf
% Y.Lu 2015/11/22 HUST/Wuhan
% input parameters
clear;clc;
K_s=0.5;            % slotting factor      
K_m=0.8;            % relative angular width of the surface mounted magnet(s)
N_sf=36;              % number of slots of full model
N_pf=6;              % number of poles of full model
N_ph=3;
p_f = N_pf/2;           % number of pole pairs
t = gcd(N_sf,N_pf/2);     %number of machine periodicity
%-------------------------单元电机中的槽极数--------------------------------
N_s=N_sf/t;
N_p=N_pf/t;
p=N_p/2;
mutual_inductance_zero =false;% The mutual inductance is not zero for default.
single_layer_wanted = 1;
% ------------winding calculations-----------------------------------------
q=N_sf/(N_ph*N_pf);   % slot per phase and per pole
yq = floor(N_s/2/p);
if (yq<1)
    yq = 1;
    mutual_inductance_zero = 0;
end


%-------可以根据自己的需要调整节距-------------------------------------------
%yq=1;节距可调整
%--------------------------------------------------------------------------


%-----------Winding feasibility check--------------------------------------
if rem(N_sf,N_ph*t) ~= 0
    fprintf('unfeasibility combination!') % N_sf/(N_ph*t) is integer and  the winding is feasible
else 
%------------Single layer feasibility and mutual inductance check----------
   if (gcd(N_sf,2)~= 2)
      single_layer_feasible = false; %! N_s is odd. SL winding is feasible only if N_s is even 
   else %! N_s is even
        if ( gcd (N_sf/t , 2 ) ~= 2 ) %!  N_s/t odd
            mutual_inductance_zero = false;
            if gcd(t,2)~=2   
                single_layer_feasible = false; %! N_s/t is odd and t is odd
            else
                single_layer_feasible = true;  %! N_s/t is odd and t is even
            end       
        else %!  N_s/t even
           if ( gcd(N_sf/(2*t), 2) ==2 )   %!  N_s / ( 2t ) even
              single_layer_feasible = true;              %!N_s/t is even and N_s/(2t) is even"; 
              if (yq==1) 
                  mutual_inductance_zero = true;  %! M = to 0 for both SL and DL winding if yq = 1
              end
           else  %!  N_s / ( 2t ) odd
                  single_layer_feasible = true;          %! N_s/t is even and N_s/(2t) is odd
               if (yq==1) 
                      mutual_inductance_zero = true;  %! M = 0 if yq = 1       
               end
               if (single_layer_wanted) 
                      mutual_inductance_zero = false;  %! M = 0 only for double layer winding
               end
           end
        end
    end                 
%end
%slot_star(); % slot_star is call
alpha_se = 2 * pi / N_s * p;         %!< Slot electrical angle*/
alpha_ph = 2 * pi / N_s *1;          %!< Angle between two adiacent phasor (or star spoke)
epsilon = 0;
if floor(N_s/(6*1)) == (N_s)/(6*1)           %!< shift of the angle if there is overlapping between spokes and the sector border
    epsilon = - alpha_ph /4;
end  
 %! First an array with all the angle and the sequence of label is create.
 %  Then the angle are sorted in order to achieve the correct sequence of spoke labels.
 %  At the same time, also the array of labels is sorted.
for i = 1:N_s         %!< Array creation
    angle(i) = alpha_se*i+epsilon;
    star(i)=i;
end
for i = 1:N_s          %!< Makes all the angles in the rang [0 2*pi]
    while (angle(i)>=2*pi)
        angle(i) = angle(i) - 2 * pi;
    end
    if (abs(angle(i)-2*pi)<0.001) 
        angle(i)=0; 
    end %!< makes 0 the angle ~ 2pi  
end
 %double swap;
 % int swap_int;
  for i = 1:N_s         %!< Sorting of the arrays
    for j = i:N_s
      if angle(j)< angle(i)
        swap     = angle(i);
        angle(i) = angle(j);
        angle(j) = swap;
        swap_int = star(i);
        star(i)  = star(j);
        star(j)  = swap_int;
      end
    end
  end
% slot_star
%! slot matrix initialization
for i = 1:N_s
    mat_A(i) = 0; %both the layer 
    mat_B(i) = 0;
    mat_C(i) = 0;
end
%! Subdivision of the star in the six sector and coil arrays population
cos_30 = cos(pi/6);  
length_A=1; length_B=1; length_C=1;  
if single_layer_wanted && single_layer_feasible %! yq must be odd
    if (gcd(yq,2)==2) 
        yq=yq-1; 
    end
    for i = 1: N_s 
      if gcd(star(i),2)==2
         continue;% skip all the even spoke
      end
      star(i) = star(i)+2;
      if star(i)>N_s star(i)=star(i)-N_s; end
         index = star(i)+yq;
      if (index>N_s) index=index-N_s;end
      
      if (cos(angle(i)) > cos_30)  
          coils_A(length_A,:) = [star(i),index,1,1,1,angle(i)];length_A=length_A+1; 
      end % A+  left side
      if (cos(angle(i)) < -cos_30) 
          coils_A(length_A,:) = [star(i),index,1,1,-1,angle(i)];length_A=length_A+1;
      end % A-
      if ( (cos(angle(i)) < 0) && (sin(angle(i)) > 0.5) ) 
          coils_B(length_B,:) = [star(i),index,1,1,1,angle(i)];length_B=length_B+1;
      end % B+  right side
      if ( (cos(angle(i)) > 0) && (sin(angle(i)) <-0.5) ) 
          coils_B(length_B,:) = [star(i),index,1,1,-1,angle(i)];length_B=length_B+1;
      end% B-
      if ( (cos(angle(i)) < 0) && (sin(angle(i)) <-0.5) ) 
          coils_C(length_C,:) = [star(i),index,1,1,1,angle(i)];length_C=length_C+1;
      end % C+
      if ( (cos(angle(i)) > 0) && (sin(angle(i)) > 0.5) ) 
          coils_C(length_C,:) = [star(i),index,1,1,-1,angle(i)];length_C=length_C+1;
      end% C-
    end     
else%double layer winding
    for i= 1: N_s
        star(i) = star(i)+2;
      if star(i)>N_s star(i)=star(i)-N_s; end
         index = star(i)+yq;
      if (index>N_s) index=index-N_s;end      
      if (cos(angle(i)) > cos_30)  
          coils_A(length_A,:) = [star(i),index,1,1,1,angle(i)];length_A=length_A+1; 
      end % A+  left side
      if (cos(angle(i)) < -cos_30) 
          coils_A(length_A,:) = [star(i),index,1,1,-1,angle(i)];length_A=length_A+1;
      end % A-
      if ( (cos(angle(i)) < 0) && (sin(angle(i)) > 0.5) ) 
          coils_B(length_B,:) = [star(i),index,1,1,1,angle(i)];length_B=length_B+1;
      end % B+  right side
      if ( (cos(angle(i)) > 0) && (sin(angle(i)) <-0.5) ) 
          coils_B(length_B,:) = [star(i),index,1,1,-1,angle(i)];length_B=length_B+1;
      end% B-
      if ( (cos(angle(i)) < 0) && (sin(angle(i)) <-0.5) ) 
          coils_C(length_C,:) = [star(i),index,1,1,1,angle(i)];length_C=length_C+1;
      end % C+
      if ( (cos(angle(i)) > 0) && (sin(angle(i)) > 0.5) ) 
          coils_C(length_C,:) = [star(i),index,1,1,-1,angle(i)];length_C=length_C+1;
      end% C-
    end  
end
 
%! Slot matrix computation %mark: slot color and text
 if single_layer_wanted && single_layer_feasible %single layer
   mark_slot = zeros(1,N_s); 
   %! Phase A  
  for i=1:(length_A-1)
    if coils_A(i,5) == 1 %sec=1
        mat_A(coils_A(i,1)) = mat_A((coils_A(i,1))) + 0.5 ;mark_slot(coils_A(i,1)) = 1; % coils_A[i].nc;
        mat_A(coils_A(i,2)) = mat_A((coils_A(i,2))) - 0.5 ;mark_slot(coils_A(i,2)) = 2;
    else
       mat_A(coils_A(i,1)) = mat_A((coils_A(i,1))) - 0.5 ; mark_slot(coils_A(i,1)) = 2;% coils_A[i].nc;
       mat_A(coils_A(i,2)) = mat_A((coils_A(i,2))) + 0.5 ; mark_slot(coils_A(i,2)) = 1;
    end
  end
  %! Phase B
  for i=1:(length_B-1)
    if coils_B(i,5) == 1 %sec=1
        mat_B(coils_B(i,1)) = mat_B(coils_B(i,1)) + 0.5 ; mark_slot(coils_B(i,1)) = 3; % coils_B[i].nc;
        mat_B(coils_B(i,2)) = mat_B(coils_B(i,2)) - 0.5 ; mark_slot(coils_B(i,2)) = 4;
    else
       mat_B(coils_B(i,1)) = mat_B(coils_B(i,1)) - 0.5 ;  mark_slot(coils_B(i,1)) = 4;% coils_B[i].nc;
       mat_B(coils_B(i,2)) = mat_B(coils_B(i,2)) + 0.5 ;  mark_slot(coils_B(i,2)) = 3;
    end
  end
  %! Phase C
  for i=1:(length_C-1)
    if coils_C(i,5) == 1 %sec=1
        mat_C(coils_C(i,1)) = mat_C(coils_C(i,1)) + 0.5 ; mark_slot(coils_C(i,1)) = 5;% coils_C[i].nc;
        mat_C(coils_C(i,2)) = mat_C(coils_C(i,2)) - 0.5 ; mark_slot(coils_C(i,2)) = 6;
    else
       mat_C(coils_C(i,1)) = mat_C(coils_C(i,1)) - 0.5 ;  mark_slot(coils_C(i,1)) = 6;% coils_C[i].nc;
       mat_C(coils_C(i,2)) = mat_C(coils_C(i,2)) + 0.5 ;  mark_slot(coils_C(i,2)) = 5;
    end
 end
 else %double layer
   mark_slot = ones(1,N_s*2); 
   for i=1:(length_A-1)
    if coils_A(i,5) == 1 %sec=1
        mat_A(coils_A(i,1)) = mat_A((coils_A(i,1))) + 0.5 ;mark_slot(coils_A(i,1)*2) = 2; % coils_A[i].nc;
        mat_A(coils_A(i,2)) = mat_A((coils_A(i,2))) - 0.5 ;mark_slot(coils_A(i,1)*2-1) = 1;
    else
       mat_A(coils_A(i,1)) = mat_A((coils_A(i,1))) - 0.5 ; mark_slot(coils_A(i,1)*2-1) = 2;% coils_A[i].nc;
       mat_A(coils_A(i,2)) = mat_A((coils_A(i,2))) + 0.5 ; mark_slot(coils_A(i,1)*2) = 1;
    end
  end  
    %! Phase B
  for i=1:(length_B-1)
    if coils_B(i,5) == 1 %sec=1
        mat_B(coils_B(i,1)) = mat_B(coils_B(i,1)) + 0.5 ; mark_slot(coils_B(i,1)*2) = 4; % coils_B[i].nc;
        mat_B(coils_B(i,2)) = mat_B(coils_B(i,2)) - 0.5 ; mark_slot(coils_B(i,1)*2-1) = 3;
    else
       mat_B(coils_B(i,1)) = mat_B(coils_B(i,1)) - 0.5 ;  mark_slot(coils_B(i,1)*2-1) = 4;% coils_B[i].nc;
       mat_B(coils_B(i,2)) = mat_B(coils_B(i,2)) + 0.5 ;  mark_slot(coils_B(i,1)*2) = 3;
    end
  end
  %! Phase C
  for i=1:(length_C-1)
    if coils_C(i,5) == 1 %sec=1
        mat_C(coils_C(i,1)) = mat_C(coils_C(i,1)) + 0.5 ; mark_slot(coils_C(i,1)*2) = 6;% coils_C[i].nc;
        mat_C(coils_C(i,2)) = mat_C(coils_C(i,2)) - 0.5 ; mark_slot(coils_C(i,1)*2-1) = 5;
    else
       mat_C(coils_C(i,1)) = mat_C(coils_C(i,1)) - 0.5 ;  mark_slot(coils_C(i,1)*2-1) = 6;% coils_C[i].nc;
       mat_C(coils_C(i,2)) = mat_C(coils_C(i,2)) + 0.5 ;  mark_slot(coils_C(i,1)*2) = 5;
    end
  end
end
mat = [mat_A;mat_B;mat_C];
% -------------------------------------------------- geometric modelling
Rm=[0.3 0.5 0.58 0.6 0.8 1]; % main radii of machine p.u.
% vizualisation settings
PoS=[0.02 0.04 0.36 0.92; 0.4 0.67 0.28 0.32; 0.4 0.34 0.28 0.32; 0.4 0.01 0.28 0.32; 0.7 0.31 0.28 0.68; 0.7 0.01 0.28 0.28];
figure(1); clf; 
scrsz = get(0,'ScreenSize'); 
set(gcf,'Position',[10 200 scrsz(3)-20 scrsz(4)/2],'Color',[1 1 1],'PaperPositionMode','auto')

ColD(1,:)=[1 0 0]; ColD(2,:)=ColD(1,:)*0.1;   % Phase A Color
ColD(3,:)=[0 1 0]; ColD(4,:)=ColD(2,:)*0.5;   % Phase B Color
ColD(5,:)=[0 0 1]; ColD(6,:)=ColD(3,:)*0.5;   % Phase C color
text_slot={'A+';'A-';'B+';'B-';'C+';'C-'};
% -------------------------------------------------- emf
wdist=[];
LinS=[':'; ':'; ':'; '-'; '-'; '-'; '-'; ];
subplot('Position',PoS(5,:)); hold on
rad=1;
%-------------------------------plot slot star-------------------------------
% for ki=1:N_s   
%     h=polar(angle(ki)*[1 1]+pi,[0 1.3]); set(h,'Color',ColD(4,:),'LineWidth',1,'LineStyle',LinS(4,:))
%     h=text(rad*cos(angle(ki)+pi),rad*sin(angle(ki)+pi),num2str(star(ki)));
%     set(h,'BackgroundColor',[1 1 .5],'Edgecolor',[.7 .7 .7],'FontSize',6,'HorizontalAlignment','center');
% end
axis([-1 1 -1 1]*1.25)
axis off, axis equal, axis tight
% -------------------------------------------------- geometry of motor parts
% subplot('Position',PoS(1,:)); hold on
% stator core
rs=[]; os=[];
if K_s==1; 
    [xs,ys]=pol2cart([0:pi/24:2*pi 2*pi:-pi/24:0],[Rm(5)*ones(1,49) Rm(6)*ones(1,49)]);
elseif K_s==0; 
    [xs,ys]=pol2cart([0:pi/24:2*pi 2*pi:-pi/24:0],[Rm(4)*ones(1,49) Rm(6)*ones(1,49)]);
else
	for ki=1:N_s
        rs=[rs Rm([5*ones(1,7) 4*ones(1,13) 5*ones(1,7)])];
        os=[os pi/N_s+[0:K_s/6:K_s K_s:(1-K_s)/6:(2-K_s) (2-K_s):K_s/6:2]*pi/N_s+(ki-1)*2*pi/N_s];
	end
	rs=[rs Rm(6)*ones(1,49)];
	os=[os (pi/N_s+2*pi):-pi/24:pi/N_s];
	[xs,ys]=pol2cart(os,rs);
end
fill(xs,ys,[0.8 0.8 0.9]-0.1,'EdgeColor',[0.8 0.8 0.9]-0.1)
%slot numeration
for ki=1:N_s
    h=text(mean(Rm(5:6))*cos(pi/N_s+(ki-1)*2*pi/N_s),mean(Rm(5:6))*sin(pi/N_s+(ki-1)*2*pi/N_s),int2str(ki));
    set(h,'BackgroundColor',[1 1 .5],'Edgecolor',[.7 .7 .7],'FontSize',8,'HorizontalAlignment','center');
end
% rotor core
[xs,ys]=pol2cart([0:pi/12:2*pi 2*pi:-pi/12:0],[Rm(1)*ones(1,25) Rm(2)*ones(1,25)]);
fill(xs,ys,[0.8 0.8 0.9]-0.1,'EdgeColor',[0.8 0.8 0.9]-0.1)
% magnets
for ki=1:N_p
    [xs,ys]=pol2cart(pi/N_p-K_m*pi/N_p+[0:pi/12:2*pi 2*pi:-pi/12:0]*K_m/N_p+(ki-1)*2*pi/N_p, ...
        [Rm(2)*ones(1,25) Rm(3)*ones(1,25)]);
    fill(xs,ys,[0.8 0.8 0.9]-0.4,'EdgeColor',[0.8 0.8 0.9]-0.4)
end
% -------------------------------------------------- geometry of windings
for ki=1:N_s
	if single_layer_wanted && single_layer_feasible %single layer
        [xs,ys]=pol2cart(pi/N_s-K_s*pi/N_s+[0:pi/12:2*pi 2*pi:-pi/12:0]*K_s/N_s+(ki-1)*2*pi/N_s, ...
        [Rm(4)*ones(1,25) Rm(5)*ones(1,25)]);
        fill(xs,ys,ColD(mark_slot(ki),:),'EdgeColor',ColD(mark_slot(ki),:));
        h=text(mean(Rm(4:5))*cos((pi/N_s-K_s*pi/N_s)+pi*K_s/N_s+(ki-1)*2*pi/N_s),...
            mean(Rm(4:5))*sin((pi/N_s-K_s*pi/N_s)+pi*K_s/N_s+(ki-1)*2*pi/N_s),text_slot(mark_slot(ki)));
        set(h,'BackgroundColor',[1 1 1],'Edgecolor',[.7 .7 .7],'FontSize',8,'HorizontalAlignment','center');  
	else
        [xs,ys]=pol2cart(-(pi/N_s-K_s*pi/N_s)-[0:pi/24:pi pi:-pi/24:0]*K_s/N_s+(ki-1)*2*pi/N_s, ...
        [Rm(4)*ones(1,25) Rm(5)*ones(1,25)]);
        fill(xs,ys,ColD(mark_slot(2*ki-1),:),'EdgeColor',ColD(mark_slot(2*ki-1),:));
        h=text(mean(Rm(4:5))*cos(-(pi/N_s-K_s*pi/N_s)-pi*K_s/N_s/2+(ki-1)*2*pi/N_s),...
            mean(Rm(4:5))*sin(-(pi/N_s-K_s*pi/N_s)-pi*K_s/N_s/2+(ki-1)*2*pi/N_s),text_slot(mark_slot(2*ki-1)));
        set(h,'BackgroundColor',[1 1 1],'Edgecolor',[.7 .7 .7],'FontSize',8,'HorizontalAlignment','center');
        [xs,ys]=pol2cart(pi/N_s-K_s*pi/N_s+[0:pi/24:pi pi:-pi/24:0]*K_s/N_s+(ki-1)*2*pi/N_s, ...
        [Rm(4)*ones(1,25) Rm(5)*ones(1,25)]);
        fill(xs,ys,ColD(mark_slot(2*ki),:),'EdgeColor',ColD(mark_slot(2*ki),:));
        h=text(mean(Rm(4:5))*cos((pi/N_s-K_s*pi/N_s)+pi*K_s/N_s/2+(ki-1)*2*pi/N_s),...
            mean(Rm(4:5))*sin((pi/N_s-K_s*pi/N_s)+pi*K_s/N_s/2+(ki-1)*2*pi/N_s),text_slot(mark_slot(2*ki)));
        set(h,'BackgroundColor',[1 1 1],'Edgecolor',[.7 .7 .7],'FontSize',8,'HorizontalAlignment','center');        
	end
end
axis equal; axis off
end