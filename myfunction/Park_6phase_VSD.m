%% 恒功率
function [Trans6pVSD]=Park_6phase_VSD(e)
%% VSD解耦变换矩阵
% Tdqz=[cos(e) cos(e-2*pi/3) cos(e+2*pi/3) cos(e-pi/6) cos(e-5*pi/6) cos(e+pi/2);
%     -sin(e) -sin(e-2*pi/3) -sin(e+2*pi/3) -sin(e-pi/6) -sin(e-5*pi/6) -sin(e+pi/2);
%     1 -1/2 -1/2 -sqrt(3)/2 sqrt(3)/2 0;
%     0 -sqrt(3)/2 sqrt(3)/2 1/2 1/2 -1;
%     1 1 1 0 0 0;
%     0 0 0 1 1 1]*(1/sqrt(3));%(这种变换只是将基波平面变换为直流量，谐波平面仍为交流量）
% ---------自然坐标系——静止坐标系----------------
Trans6pVSD.T6s=1/sqrt(3)*[1 -1/2 -1/2 sqrt(3)/2 -sqrt(3)/2 0;
    0 sqrt(3)/2 -sqrt(3)/2 1/2 1/2 -1;
    1 -1/2 -1/2 -sqrt(3)/2 sqrt(3)/2 0;
    0 -sqrt(3)/2 sqrt(3)/2 1/2 1/2 -1;
    1 1 1 0 0 0;
    0 0 0 1 1 1];
% ---------静止坐标系——同步旋转坐标系----------------
Trans6pVSD.Ts2r=[cos(e) sin(e) 0 0 0 0;
    -sin(e) cos(e) 0 0 0 0;
    0 0 -cos(e) sin(e) 0 0;
    0 0 sin(e) cos(e) 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];%（这种变换将基波平面和谐波平面均变换为直流量）
% ---------自然坐标系——同步旋转坐标系----------------
Trans6pVSD.Tdqz=Trans6pVSD.Ts2r*Trans6pVSD.T6s;
end