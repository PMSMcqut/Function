clear;clc;
%% Unit
mm=1e-3;
mu0=4*pi*10^-7;% permeability in vacuum
Br=load('Br-rectangular.txt');
Bt=load('Bt-rectangular.txt');
p=3;
Qs=36;
Nt=gcd(Qs,p);% number of motor units
rpm=1000;% speed 
f0=p*rpm/60;% fundamental electrical frequncy
T=1/f0;% electrical period
TimeStep=200;
AngleStep=360;
Time=0:T/TimeStep:T;% Time for FFT
Space=0:360/AngleStep:360;
TimeDomain.Time=Time(1:TimeStep);
TimeDomain.Space=Space(1:AngleStep);% Space angle for FFT
%% Harmonic analysis
% ====================== flux density ==============================
Br=Br(1:AngleStep,1:TimeStep);
Bt=Bt(1:AngleStep,1:TimeStep);
[FourierBr]=FFT_fun(Br,TimeDomain.Time,TimeDomain.Space,f0,5000,'2D','fun');
[FourierBt]=FFT_fun(Bt,TimeDomain.Time,TimeDomain.Space,f0,5000,'2D','fun');
% ====================== force density =============================
F.rad=(Br.^2-Bt.^2)/(2*mu0);
F.tan=Br.*Bt/mu0;
[FourierFr]=FFT_fun(F.rad,TimeDomain.Time,TimeDomain.Space,f0,5000,'2D','fun');
[FourierFt]=FFT_fun(F.tan,TimeDomain.Time,TimeDomain.Space,f0,5000,'2D','fun');
%% ================= FigPlot ========================================
% figure(1)
% FluxPlot(Br,Bt,FourierBr,FourierBt,TimeDomain.Time,TimeDomain.Space);
figure(1)
ForcePlot2(F,FourierFr,FourierFt,TimeDomain.Time,TimeDomain.Space);
% figure(3)
% Harmonic2D(Br,Bt,F,Time,TimeStep,Space,AngleStep,f0);