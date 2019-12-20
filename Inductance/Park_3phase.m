%% 恒功率
function [Trans3p]= Park_3phase (e)
% ---------三相自然坐标系——两相静止坐标系----------------
Trans3p.T3s2s=sqrt(2/3)*[1 -1/2 -1/2;
    0 sqrt(3)/2 -sqrt(3)/2;
    sqrt(2)/2 sqrt(2)/2 sqrt(2)/2];
    % ---------两相静止坐标系——同步旋转坐标系----------------
    Trans3p.T2s2r=[cos(e) sin(e) 0;
        -sin(e) cos(e) 0;
        0 0 1];
    % ---------三相自然坐标系——同步旋转坐标系----------------
    Trans3p.T3s2r=Trans3p.T2s2r*Trans3p.T3s2s;
end