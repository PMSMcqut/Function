%%%%%%%%%%%%%
%本程序旨在与分析IPM(包括SPM)电机的过载性能，主要由两部分组成：1.额定MTPA FW，MTPV；2过载情况下MTPA FW，MTPV；
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输入五个输入量，psim ld，lq，i，u都是标幺值  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
psim=input('please input psim pu:');
Ld=input('please input Ld pu:');
Lq=input('please input Lq pu:');
I1=input('please input I pu:');
U=input('please input U pu:');
%psim=0.6;
%Ld=0.3;
%Lq=0.9;
%I1=1;
%U=1;
speed=0;
gamma_max=asin(0.25*(psim/(Ld-Lq)/I1+sqrt((psim/(Ld-Lq)/I1)^2+8)));
speed_corner=U/sqrt((psim-Ld*I1*sin(gamma_max))^2+(Lq*I1*cos(gamma_max))^2);
speed_corner_SPM=U/sqrt((Ld*I1)^2+(psim)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=0;
if Ld==Lq
    if psim/Ld>=I1-0.0001
        while speed<speed_corner_SPM
           n=n+1;
           gamma(n)=0;
           Te(n)=psim*I1*cos(gamma(n));
           rev(n)=speed;
           speed=speed+0.05;
        end
        while speed>=speed_corner_SPM
            n=n+1;
            rev(n)=speed;
            gamma(n)=asin(((Ld*I1)^2+psim^2-(U/rev(n))^2)/(2*Ld*I1*psim));
            Te(n)=psim*I1*cos(gamma(n));

            speed=speed+0.05;
            if rev(n)>5
                break
            end
        end

    else 
         while speed<speed_corner_SPM
           n=n+1;
           gamma(n)=0;
           Te(n)=psim*I1*cos(gamma(n));
           rev(n)=speed;
           speed=speed+0.05;
        end
        while speed>=speed_corner_SPM
            n=n+1;
            rev(n)=speed;
            gamma(n)=asin(((Ld*I1)^2+psim^2-(U/rev(n))^2)/(2*Ld*I1*psim));
            Te(n)=psim*I1*cos(gamma(n));
            speed=speed+0.05;
            if I1*sin(gamma(n))>=psim/Ld%%进入MTPV zone
                while rev(n)<5
                n=n+1;
                rev(n)=speed;
                Te(n)=psim*U/rev(n)/Ld;
                speed=speed+0.05;

                end
            if rev(n)>5
                break
            end
                
            end

        end
    end
plot(rev,Te.*rev,'g^-');
grid on
hold on
clear
else
    if psim/Ld>=I1-0.001
    while speed<speed_corner 
        n=n+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MTPA zone gamma 是弧度
%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma(n)=asin(0.25*(psim/(Ld-Lq)/I1+sqrt((psim/(Ld-Lq)/I1)^2+8)));
Te(n)=psim*I1*cos(gamma(n))-0.5*(Ld-Lq)*I1^2*sin(2*gamma(n)); 
rev(n)=speed;
speed=speed+0.05;
end
 
while speed>=speed_corner
n=n+1;
k=(2*Ld*I1*psim+sqrt((2*Ld*I1*psim)^2-4*((Ld*I1)^2-(Lq*I1)^2)*(psim^2+(Lq*I1)^2-(U/speed)^2)))/(2*((Ld*I1)^2-(Lq*I1)^2));
if k>0    
gamma(n)=asin(k);
end
if k<0
 gamma(n)=asin((2*Ld*I1*psim-sqrt((2*Ld*I1*psim)^2-4*((Ld*I1)^2-(Lq*I1)^2)*(psim^2+(Lq*I1)^2-(U/speed)^2)))/(2*((Ld*I1)^2-(Lq*I1)^2)));
end
Te(n)=psim*I1*cos(gamma(n))-0.5*(Ld-Lq)*I1^2*sin(2*gamma(n)); 
rev(n)=speed;
speed=speed+0.05;
if rev(n)>5
    break
end
  end
plot(rev,Te.*rev,'g^-');
grid on
hold on
clear
else
   while speed<speed_corner 
        n=n+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MTPA zone gamma 是弧度
%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma(n)=asin(0.25*(psim/(Ld-Lq)/I1+sqrt((psim/(Ld-Lq)/I1)^2+8)));
Te(n)=psim*I1*cos(gamma(n))-0.5*(Ld-Lq)*I1^2*sin(2*gamma(n)); 
rev(n)=speed;
speed=speed+0.05;
end
 
while speed>=speed_corner 
n=n+1;
k=(2*Ld*I1*psim+sqrt((2*Ld*I1*psim)^2-4*((Ld*I1)^2-(Lq*I1)^2)*(psim^2+(Lq*I1)^2-(U/speed)^2)))/(2*((Ld*I1)^2-(Lq*I1)^2));
if k>0    
gamma(n)=asin(k);
end
if k<0
 gamma(n)=asin((2*Ld*I1*psim-sqrt((2*Ld*I1*psim)^2-4*((Ld*I1)^2-(Lq*I1)^2)*(psim^2+(Lq*I1)^2-(U/speed)^2)))/(2*((Ld*I1)^2-(Lq*I1)^2)));
end
    speed=speed+0.01;
    rev(n)=speed;
    Te(n)=psim*I1*cos(gamma(n))-0.5*(Ld-Lq)*I1^2*sin(2*gamma(n)); 
    z=rev(n)*(psim/U)/(1-Ld/Lq);
    s=acos((z-sqrt(z^2+8))/4);
    
   %speed
   %sqrt(((U*cos(s)-speed*psim)/speed/Ld)^2+(U*sin(s)/speed/Lq))
   % Id=I1*sin(gamma(n));
   %Id=(U*cos(s)-speed*psim)/speed/Ld;
   %Iq=I1*cos(gamma(n))
   %Iq=U*sin(s)/speed/Lq
   %Te(n)=psim*U*sin(s)/speed/Lq+(Ld-Lq)*(U*cos(s)-speed*psim)/speed/Ld*U*sin(s)/speed/Lq
   if I1>=sqrt(((U*cos(s)-speed*psim)/speed/Ld)^2+(U*sin(s)/speed/Lq)^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%进入MTPV zone gamma 是弧度,s为内功率因数角
%%%%%%%%%%%%%%%%%%%%%%%%%%%
while rev(n)<5;%%%%%以最大转速作为MTPV的终结篇
       Te(n)=psim*U*sin(s)/speed/Lq+(Ld-Lq)*(U*cos(s)-speed*psim)/speed/Ld*U*sin(s)/speed/Lq;
      %Te(n)=1/rev(n)*(U*psim/Ld*sin(s)+1/2*U^2*(1/Lq-1/Ld)*sin(2*s));
      n=n+1;
      Te(n)=psim*U*sin(s)/speed/Lq+(Ld-Lq)*(U*cos(s)-speed*psim)/speed/Ld*U*sin(s)/speed/Lq;
      %Te(n)=1/rev(n-1)*(U*psim/Ld*sin(s)+1/2*U^2*(1/Lq-1/Ld)*sin(2*s));
      speed=speed+0.05;
      rev(n)=speed;
      z=rev(n)*(psim/U)/(1-Ld/Lq);
      s=acos((z-sqrt(z^2+8))/4);
     

end
if rev(n)>5
break
end
    
%if gamma(n)>1.5
    %break
%end

  end
end
plot(rev,Te.*rev,'r*-');
grid on
hold on
clear
end

end

