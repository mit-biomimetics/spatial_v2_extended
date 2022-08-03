function [p_link,R_link,X0_link] = get_link_position( model, q, link_idx)
% p_link = get_link_position( model, q, gc) computes the position "p_link"
% (in the world frame) of the origin of frame corresponding to the 
% robot's link with index "link_idx" for the given robot "model" in 
% the given joint configuration "q".
%% Initialization

X0  = cell(model.NB,1);
p0  = cell(model.NB,1);
R0  = cell(model.NB,1);
%% floating base forward kinematics

if ~strcmp(model.fb_type,'eul')
    error('get_gc_position only works with euler angle-based floating base joint')
end

[dim_fb, Xup] = fwd_kin_fb(model,q);

%% Joints Forward Kinematics

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
