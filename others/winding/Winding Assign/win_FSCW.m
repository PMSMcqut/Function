% calculation of three-phase winding distribution and mmf
% Jian Li 2011 20th Jan. PEAL/DAU 
% input parameters
function [mat,mat_A,mat_B,mat_C,k_w,single_layer_feasible,t,winding_feasibility,q]=win_FSCW(N_sf,N_pf,N_ph,single_layer_wanted)
p_f = N_pf/2;           % number of pole pairs
yq = -1;
t = gcd(N_sf,N_pf/2);     %number of machine periodicity

mutual_inductance_zero = false;% The mutual inductance is not zero for default.
%single_layer_wanted = 1;
%Winding feasibility
winding_feasibility = 1;
if rem(N_sf,N_ph*t) ~= 0
    disp('unfeasibility combination!')
    winding_feasibility = 0;
    %return;
end
% -------------------------------------------------- winding calculations
q=N_sf/(N_ph*N_pf);   % slot per phase and per pole
%alpha_ph=2*pi/(N_s/t);% angle between two spokes
%alpha_se=alpha_ph*N_p/2/t;%angle between the phasors of two adjacent slots in electrical angle
%N_s=N_sf/t;
%N_p=N_pf/t;
%p=N_p/2;
N_s=N_sf;
N_p=N_pf;
p=N_pf/2;

yq = floor(N_s/2/p);
%if q <1
%    yq = 1;
%elseif q>=1
%    yq=fix(N_s,N_p)
%end
if (yq<1)
    yq = 1;
    mutual_inductance_zero = 0;
end
%Winding feasibility check
if ( N_sf/(3*t) ~= (N_sf+0.0)/(3*t)) %! N_s/ (3t) is not integer
  t = -1;
else  %!<  N_s/ (3t) is integer and  the winding is feasible
%! Single layer feasibility and mutual inductance check
  if (gcd(N_sf,2)~= 2)
  single_layer_feasible = false; %! N_s is odd. SL winding is feasible only if N_s is even 
  else %! N_s is even
        if ( gcd (N_sf/t , 2 ) ~= 2 ) %!  N_s/t odd
            mutual_inductance_zero = false;
            if gcd(t/2,2)~=2   
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
%slot_star(); % slot_star is call
end

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


ColD(1,:)=[1 0 0]; ColD(2,:)=ColD(1,:)*0.1;   % Phase A Color
ColD(3,:)=[0 1 0]; ColD(4,:)=ColD(2,:)*0.5;   % Phase B Color
ColD(5,:)=[0 0 1]; ColD(6,:)=ColD(3,:)*0.5;   % Phase C color
text_slot={'A+';'A-';'B+';'B-';'C+';'C-'};

 for n=1:1000
    k_w(n)=kw(n,mat_A,N_s,t);
 end


