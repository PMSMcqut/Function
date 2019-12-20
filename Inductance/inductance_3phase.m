function [Ldq0]=inductance_3phase(Lph,e)
% -------------------- Fourier transfer------------------------------
Temp=[Lph(:,1),Lph(:,2)];
[Spectrum]=FFT_fun(Temp,time,1,1,1,'1D','fun');
% -------------------- Inductance Calcualtion------------------------
for i=1:length(e) % ´ú±íÊ±¿Ì
    Ls=[Lph(i,1) Lph(i,2) Lph(i,3);
        Lph(i,4) Lph(i,5) Lph(i,6);
        Lph(i,7) Lph(i,8) Lph(i,9)];
    [Trans3p]= Park_3phase (e);
    Ldq0=Trans3p.T3s2r*Ls*inv(Trans3p.T3s2r);
    Ld(i)=Ldqzi(1,1);
Lq(i)=Ldqzi(2,2);
L0(i)=Ldqzi(3,3);
end
Ldq0=[mean(Ld) 0 0;
    0 mean(Lq) 0;
    0 0 mean(L0)];
end