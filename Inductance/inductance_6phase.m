function [Spectrum,Lana,L2dq,Lvsd]=inductance_6phase(Lph,time,e)
Temp=[Lph(:,1),Lph(:,2),Lph(:,4)];
% -------------------- Fourier transfer------------------------------
[Spectrum]=FFT_fun(Temp,time,1,1,1,'1D','fun');
% -------------------- analytical method ----------------------------
Lana.Ls0=Spectrum.P.Amplitude(1,1);%自感直流分量,A-A
Lana.Ls2=-Spectrum.P.Amplitude(3,1);%自感二次分量,A-A
Lana.Ms0=-Spectrum.P.Amplitude(1,2);%一个三相内直流分量.A-B
Lana.Ms2=-Spectrum.P.Amplitude(3,2);%一个三相内二次分量,A-B
Lana.Mm0=Spectrum.P.Amplitude(1,3)/(sqrt(3)/2);%两个三相间直流分量,A-D  
Lana.Mm2=-Spectrum.P.Amplitude(3,3);%两个三相间二次分量,A-D

Lana.Ld1=Lana.Ls0+1/2*Lana.Ls2-Lana.Ms0+Lana.Ms2;
Lana.Lq1=Lana.Ls0-1/2*Lana.Ls2-Lana.Ms0-Lana.Ms2;
Lana.Md1d2=3/2*(Lana.Mm0+Lana.Mm2);
Lana.Mq1q2=3/2*(Lana.Mm0-Lana.Mm2);

Lana.Ld=Lana.Ld1+Lana.Md1d2;
Lana.Lq=Lana.Lq1+Lana.Mq1q2;
Lana.Ldz=Lana.Ld1-Lana.Md1d2;
Lana.Lqz=Lana.Lq1-Lana.Mq1q2;
% Lana.Ld=1/2*(2*Lana.Ls0 + Lana.Ls2 + 3*Lana.Mm0 + 3*Lana.Mm2 - 2*Lana.Ms0 + 2*Lana.Ms2);
% Lana.Lq=1/2*(2*Lana.Ls0 - Lana.Ls2 + 3*Lana.Mm0 - 3*Lana.Mm2 - 2*Lana.Ms0 - 2*Lana.Ms2);
% Lana.Ldz=1/2*(2*Lana.Ls0 + Lana.Ls2 - 3*Lana.Mm0 - 3*Lana.Mm2 - 2*Lana.Ms0 + 2*Lana.Ms2);
% Lana.Lqz=1/2*(2*Lana.Ls0 - Lana.Ls2 - 3*Lana.Mm0 + 3*Lana.Mm2 - 2*Lana.Ms0 - 2*Lana.Ms2);

for i=1:length(e) % 代表时刻
L11=[Lph(i,1) Lph(i,2) Lph(i,3);
    Lph(i,7) Lph(i,8) Lph(i,9);
    Lph(i,13) Lph(i,14) Lph(i,15)];
M12=[Lph(i,4) Lph(i,5) Lph(i,6);
    Lph(i,10) Lph(i,11) Lph(i,12);
    Lph(i,16) Lph(i,17) Lph(i,18)];
M21=[Lph(i,19) Lph(i,20) Lph(i,21);
    Lph(i,25) Lph(i,26) Lph(i,27);
    Lph(i,31) Lph(i,32) Lph(i,33)];
L22=[Lph(i,22) Lph(i,23) Lph(i,24);
    Lph(i,28) Lph(i,29) Lph(i,30);
    Lph(i,34) Lph(i,35) Lph(i,36)];
Ls=[L11 M12;
    M21 L22];
% --------------------- double d-q transform --------------------------
[Trans6p2dq]=Park_6phase_double_dq(e(i));
Ldq1=Trans6p2dq.T3s2r1*L11*inv(Trans6p2dq.T3s2r1);
Mdq12=Trans6p2dq.T3s2r1*M12*inv(Trans6p2dq.T3s2r2);
Ldq2=Trans6p2dq.T3s2r2*L22*inv(Trans6p2dq.T3s2r2);
Mdq21=Trans6p2dq.T3s2r2*M21*inv(Trans6p2dq.T3s2r1);
Ld1(i)=Ldq1(1,1);Lq1(i)=Ldq1(2,2);
Md12(i)=Mdq12(1,1);Mq12(i)=Mdq12(2,2);
Ld2(i)=Ldq2(1,1);Lq2(i)=Ldq2(2,2);
Md21(i)=Mdq21(1,1);Mq21=Mdq21(2,2);
% --------------------- VSD transform ------------------------------
[Trans6pVSD]=Park_6phase_VSD(e(i));
Ldqzi=Trans6pVSD.Tdqz*Ls*inv(Trans6pVSD.Tdqz);
Ld(i)=Ldqzi(1,1);
Lq(i)=Ldqzi(2,2);
Ldz(i)=Ldqzi(3,3);
Lqz(i)=Ldqzi(4,4);
Lo1(i)=Ldqzi(5,5);
Lo2(i)=Ldqzi(6,6);
end
% --------------------- double d-q transform --------------------------
L2dq.L0=[mean(Ld1) 0 mean(Md12) 0;
    0 mean(Lq1) 0 mean(Mq12);
    mean(Md21) 0 mean(Ld2) 0;
    0 mean(Mq21) 0 mean(Lq2);];
L2dq.Ldqz=Trans6p2dq.T_decouple*L2dq.L0*inv(Trans6p2dq.T_decouple);
% --------------------- VSD transform ------------------------------
Lvsd.dqz=[mean(Ld) 0 0 0 0 0;
    0 mean(Lq) 0 0 0 0;
    0 0 mean(Ldz) 0 0 0;
    0 0 0 mean(Lqz) 0 0;
    0 0 0 0 mean(Lo1) 0;
    0 0 0 0 0 mean(Lo2);];
end
