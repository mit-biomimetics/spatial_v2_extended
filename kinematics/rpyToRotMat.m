function R_body_to_world = rpyToRotMat(rpy)
% Need to define this rotation better 
R_body_to_world = rz(rpy(3))' * ry(rpy(2))' * rx(rpy(1))'; 