% Inductance calculation for PMSM
% Copyright (c) 2017-2020
% School of Electrical and Electronics Engineering
% Huazhong University of Science and Technology
% All Rights Reserved
% 
%              Programmer:Lu,Yang yanglu_hust@sina.com 
%              Version:2.0
%              Last modified:11/29/2018
%--------------------------introduction-------------------------
% The data is the flux linkage of each phase,A-A,A-B,A-C,A-D£¬A-E,A-F.....,F-F
% Lana:     analytical method
% L2dq:     double d-q method
% Lvsd:     VSD method
clear;clc;
%% data read
load('current.txt');
Iph=current(:,2:size(current,2));
time=current(:,1);
flux=load('fluxlinkage.txt');
phase=6;
p=3;
rpm=1500;
f=p*rpm/60;
e=2*pi*f*time;
%% data process
% -------------------- phase inductance calculation,Unit: mH------------------------------
for i=1:phase
    for j=1:phase
        Lph(:,j+phase*(i-1))=flux(:,j+phase*(i-1))./Iph(:,i)*1000;
    end
end
%% calculation transfer martrix
% ----------------´´½¨º¯Êý¾ä±ú ---------------------------
% Park_3p=@Park_3phase;
% Park_6pVSD=@Park_6phase_VSD;
% Park_6p2dq=@Park_6phase_double_dq;
switch phase
    case 3
        [Ldq0]=inductance_3phase(Lph,e);
    case 6
        [Spectrum,Lana,L2dq,Lvsd]=inductance_6phase(Lph,time,e);
end


