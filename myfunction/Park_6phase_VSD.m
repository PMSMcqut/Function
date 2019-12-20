%% �㹦��
function [Trans6pVSD]=Park_6phase_VSD(e)
%% VSD����任����
% Tdqz=[cos(e) cos(e-2*pi/3) cos(e+2*pi/3) cos(e-pi/6) cos(e-5*pi/6) cos(e+pi/2);
%     -sin(e) -sin(e-2*pi/3) -sin(e+2*pi/3) -sin(e-pi/6) -sin(e-5*pi/6) -sin(e+pi/2);
%     1 -1/2 -1/2 -sqrt(3)/2 sqrt(3)/2 0;
%     0 -sqrt(3)/2 sqrt(3)/2 1/2 1/2 -1;
%     1 1 1 0 0 0;
%     0 0 0 1 1 1]*(1/sqrt(3));%(���ֱ任ֻ�ǽ�����ƽ��任Ϊֱ������г��ƽ����Ϊ��������
% ---------��Ȼ����ϵ������ֹ����ϵ----------------
Trans6pVSD.T6s=1/sqrt(3)*[1 -1/2 -1/2 sqrt(3)/2 -sqrt(3)/2 0;
    0 sqrt(3)/2 -sqrt(3)/2 1/2 1/2 -1;
    1 -1/2 -1/2 -sqrt(3)/2 sqrt(3)/2 0;
    0 -sqrt(3)/2 sqrt(3)/2 1/2 1/2 -1;
    1 1 1 0 0 0;
    0 0 0 1 1 1];
% ---------��ֹ����ϵ����ͬ����ת����ϵ----------------
Trans6pVSD.Ts2r=[cos(e) sin(e) 0 0 0 0;
    -sin(e) cos(e) 0 0 0 0;
    0 0 -cos(e) sin(e) 0 0;
    0 0 sin(e) cos(e) 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];%�����ֱ任������ƽ���г��ƽ����任Ϊֱ������
% ---------��Ȼ����ϵ����ͬ����ת����ϵ----------------
Trans6pVSD.Tdqz=Trans6pVSD.Ts2r*Trans6pVSD.T6s;
end