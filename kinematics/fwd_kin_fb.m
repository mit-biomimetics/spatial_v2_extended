function [Xup, S]  = fwd_kin_fb(model,q)
% this function provides the floating base dimensions and floating base
% coordinates for a 3D or planar spatial model. The type of model is
% determined by the field calle 'type' in the model structure. The value of
% model.type must be a string 'planar' or '3D'. If no field type is
% specified, the default is a 3D model.

S = cell(model.NB,1);

if strcmp(model.fb_type, 'planar')
    rpy = [0;q(3);0];
    tr = [q(1);0;q(2)];
    S{3} = [
            0 0 0;
            1 0 0;
            0 0 0;
            0 1 0;
            0 0 0;
            0 0 1];
elseif strcmp(model.fb_type,'eul') % 3D floating base
    rpy = q(4:6);
    tr = q(1:3);
    S{6} = eye(6);
else
    error('fb type must be "planar" or "eul" ')
end


Xup = cell(model.NB, 1);
for i = 1:(model.NB_fb - 1)
    Xup{i} = zeros(6,6);
end

R_world_to_body = rpyToRotMat(rpy)';
Xup{model.NB_fb} = [R_world_to_body zeros(3,3);...
    -R_world_to_body*skew(tr) R_world_to_body];
