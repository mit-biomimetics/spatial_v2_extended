function R_body_to_world = rpyToRotMat(rpy)
% Need to define this rotation better

if (size(rpy) == [3 1])
    R_body_to_world = rz(rpy(3))' * ry(rpy(2))' * rx(rpy(1))';
elseif (size(rpy) == [1 1])
    R_body_to_world = [1 0 0; 0 0 1] * ry(rpy)' * [1 0 0; 0 0 1]';
else
    error('Invalid size for rpy')
end
