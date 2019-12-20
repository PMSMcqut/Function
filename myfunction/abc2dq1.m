
% abc -> dq general transformation for multi 3-phase machines
% same results as dq2abc.m for the case n3phase=1

function idq = abc2dq(i123, theta_i)
i123 = i123';
% matrice di trasformazione (3 -> 2)
T32_d = 2/3 * [cos(theta_i)     cos(theta_i-2*pi/3)    cos(theta_i-4*pi/3)];
T32_q = 2/3 * [-sin(theta_i)    -sin(theta_i-2*pi/3)   -sin(theta_i-4*pi/3)];
% 123 -> alpha beta

for i = 1 : size(i123 , 2)
    id(i, 1) = T32_d(i,:)*i123(:,i);
    iq(i, 1) = T32_q(i,:)*i123(:,i);
end

idq = [id iq];


