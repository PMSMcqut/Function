%% 恒功率
function [Trans6p2dq]=Park_6phase_double_dq(e)
% ---------三相自然坐标系——两相静止坐标系----------------
Trans6p2dq.T3s2s=sqrt(2/3)*[1 -1/2 -1/2;
    0 sqrt(3)/2 -sqrt(3)/2;
    sqrt(2)/2 sqrt(2)/2 sqrt(2)/2];
%% 第一套三相绕组变换矩阵
    % ---------两相静止坐标系——同步旋转坐标系----------------
    Trans6p2dq.T2s2r1=[cos(e) sin(e) 0;
        -sin(e) cos(e) 0;
        0 0 1];
    % ---------三相自然坐标系——同步旋转坐标系----------------
    Trans6p2dq.T3s2r1=Trans6p2dq.T2s2r1*Trans6p2dq.T3s2s;
    
%% 第二套三相绕组变换矩阵
    % ---------两相静止坐标系——同步旋转坐标系----------------
    Trans6p2dq.T2s2r2=[cos(e-pi/6) sin(e-pi/6) 0;
        -sin(e-pi/6) cos(e-pi/6) 0;
        0 0 1];
    % ---------三相自然坐标系——同步旋转坐标系----------------
    Trans6p2dq.T3s2r2=Trans6p2dq.T2s2r2*Trans6p2dq.T3s2s;
%% 解耦矩阵
    Trans6p2dq.T_decouple=1/2*[1 0 1 0;
        0 1 0 1;
        -1 0 1 0;
        0 1 0 -1];
end