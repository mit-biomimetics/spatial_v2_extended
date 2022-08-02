function [dim_fb,tr,rpy]  = get_fb_dim(model,q)
% this function provides the floating base dimensions and floating base
% coordinates for a 3D or planar spatial model. The type of model is
% determined by the field calle 'type' in the model structure. The value of 
% model.type must be a string 'planar' or '3D'. If no field type is
% specified, the default is a 3D model.
if isfield(model,'type')
    if strcmp(model.type,'planar')
        dim_fb = 3;
        rpy = [0;q(3);0];
        tr = [q(1);0;q(2)];
    elseif strcmp(model.type,'3D') % 3D floating base
        dim_fb = 6;
        rpy = q(4:6);
        tr = q(1:3);
    else 
        error('type must be "planar" or "3D" ')
    end
else % 3D floating base
    dim_fb = 6;
    rpy = q(4:6);
    tr = q(1:3);
end