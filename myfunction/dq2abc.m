% dq -> abc general transformation for multi 3-phase machines

function iabc = dq2abc(id, iq, theta_i)

% Inverse transformation matrix (2 -> 3)
T23 = [cos(theta_i)            -sin(theta_i);
       cos(theta_i - 2*pi/3)     -sin(theta_i - 2*pi/3);
       cos(theta_i + 2*pi/3)     -sin(theta_i + 2*pi/3)];
   
% idq = (id + 1j* iq);
% d q -> 123
iabc = T23 * [real(idq);imag(idq)];
