%input
nmax=11000;%rpm
Tmax=310; %N
Pmax=150; %kw
TN=143;%N
PN=60; %kw
%calc
nc=9549*Pmax/Tmax; %rpm
nN=9549*PN/TN; %rpm
%dim curve
Tcont=ones(1,1000);
Pcont=ones(1,1000);
Tpk=ones(1,1000);
Ppk=ones(1,1000);
%calc continues courve
for k=1:1:1000
    n(k)=nmax/1000*k;
    if n(k)<=nN
        Tcont(k)=TN;
    else
        Tcont(k)=9549*PN/n(k);
    end
    Pcont(k)=Tcont(k)*n(k)/9549;
end
%calc peak courve
for k=1:1:1000
    if n(k)<=nc
        Tpk(k)=Tmax;
    else
        Tpk(k)=9549*Pmax/n(k);
    end
    Ppk(k)=Tpk(k)*n(k)/9549;
end
%plot curve
figure(1) 
%plot(n,Tcont,'--',n,Pcont,'--',n,Tpk,'-',n,Ppk,'-');
yyaxis left;%******  
plot(n,Tcont,'b--',n,Tpk,'b-');
xlabel('转速/rpm');  
ylabel('转矩/Nm');  
yyaxis right;%******  
plot(n,Pcont,'r--',n,Ppk,'r-');
ylabel('功率/kW');  

title('电机外特性');
legend('连续转矩','连续功率','峰值转矩','峰值功率')        
    
