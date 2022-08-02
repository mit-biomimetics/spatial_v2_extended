function [p_link,R_link,X0_link] = get_link_position( model, q, link_idx)
% p_link = get_link_position( model, q, gc) computes the position "p_link"
% (in the world frame) of the origin of frame corresponding to the 
% robot's link with index "link_idx" for the given robot "model" in 
% the given joint configuration "q".
%% Get floating base dimensions (for 3D model support)
if ~strcmp(model.fb_type,'eul')
    error('get_gc_position only works with euler angle-based floating base joint')
end
[dim_fb, tr, rpy] = get_fb_dim(model,q);

%%
Xup = cell(model.NB,1);
X0  = cell(model.NB,1);
p0  = cell(model.NB,1);
R0  = cell(model.NB,1);

%% Forward Kinematics
R_world_to_body = rpyToRotMat(rpy)';

for i = 1:(dim_fb - 1)
    Xup{i} = zeros(6,6);
end

Xup{dim_fb} = [R_world_to_body zeros(3,3);...
         -R_world_to_body*skew(tr) R_world_to_body];

for i = (dim_fb + 1):model.NB
    [ XJ, ~ ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
end
X0{dim_fb} = Xup{dim_fb};

for i = (dim_fb+1):model.NB
    X0{i} = Xup{i} * X0{model.parent(i)};                      % propagate xform from origin
    [R0{i}, p0{i}] = plux(X0{i});                                % rotation from origin, translation from origin
end

for i = 1:length(link_idx)
    indx = (3*i-2):(3*i);
    p_link(indx) = p0{link_idx(i)};
    R_link(:,:,i) = R0{link_idx(i)};
end

X0_link = {X0{link_idx}};
