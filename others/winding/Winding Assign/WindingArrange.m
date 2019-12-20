%% calculating the slot number of phase shift
function [Winding,Slot]=WindingArrange(Qs,p,yq,m,n3ph,Beta0,alphaS)
for DeltaS=1:Qs
    if mod(DeltaS*p*360/Qs+alphaS,360)==0||mod(DeltaS*p*360/Qs-alphaS,360)==0
        break
    end
end
%% begin to caculate the star of slots
Qs0=Qs/n3ph;
m0=360/Beta0;
g=gcd(Qs0,6*p); % caculate the gratest common division of slots and 2*m*p
N=Qs0/g;D=6*p/g;
if Beta0==60
    K=0;
else
    K=N;
end
J=[1:Qs0,-Qs0:-1];% the number of 3-phase slots
%% ------caculate the position of each slot
R=D*(abs(J)-1)+1+3*N*((1-sign(J))/2)-6*N;
for i=1:1:2*Qs0
    while R(i)>6*N
     R(i)=R(i)-6*N;
    end
    while R(i)<=0
        R(i)=R(i)+6*N;
    end
end
j=0;
%% distribution the slots to A B C phase
for i=1:1:2*Qs0
   if R(i)<=(N+K*(1+sign(J(i)))/2)&&R(i)>(K*(1-sign(J(i)))/2)
      j=j+1; A(j)=J(i);
   end
   if R(i)>(2*N+K*(1-sign(J(i)))/2)&&R(i)<=(3*N+K*(1+sign(J(i)))/2)
      j=j+1; B(j)=J(i);
   end
   if R(i)>(4*N+K*(1-sign(J(i)))/2)&&R(i)<=(5*N+K*(1+sign(J(i)))/2)
      j=j+1; C(j)=J(i);
   end
   
end
A(A==0)=[];
B(B==0)=[];
C(C==0)=[];
Winding_temp=[A;B;C];
for ii=1:n3ph
    WindingTop(ii,:)=sort([A,B,C],'ComparisonMethod','abs');% 将ABC三相按照槽号排序
    for i=1:3
         [lia,locb]=ismember(Winding_temp(i,:),WindingTop(ii,:));%ABC三相绕组上层边所处的位置
         WindingTop(ii,locb)=i*sign(WindingTop(ii,locb));%用1，2，3代表A，B，C三相绕组,4,5,6代表D,E,F绕组
    end
    WindingTop(ii,:)=WindingTop(ii,:)+sign(WindingTop(ii,:))*3*(ii-1);
end
if n3ph==1
    Winding.top=WindingTop;
else
for j=1:Qs0
    Winding.top(2*j-1)=WindingTop(1,j);
    if 2*j-1+DeltaS<=Qs
        Winding.top(2*j-1+DeltaS)=WindingTop(2,j);
    else
        Winding.top(abs(2*j-1+DeltaS-Qs))=WindingTop(2,j);
    end
end
end
for j=1:Qs
      if j+yq<=Qs
        Winding.Bottom(j+yq)=-Winding.top(j);
      else
        Winding.Bottom(abs(j+yq-Qs))=-Winding.top(j);
      end
end
Winding.all=[Winding.top;Winding.Bottom];
for iii=1:m
    Slot(iii,:)=find(abs(Winding.all(1,:))==iii).*sign(Winding.all(1,(find(abs(Winding.all(1,:))==iii))));
end

end