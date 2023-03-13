function [p_link,R_link,X0_link] = get_link_position( model, q, link_idx)
% p_link = get_link_position( model, q, gc) computes the position "p_link"
% (in the world frame) of the origin of frame corresponding to the 
% robot's link with index "link_idx" for the given robot "model" in 
% the given joint configuration "q".
if ~strcmp(model.fb_type,'eul')
    error('get_gc_position only works with euler angle-based floating base joint')
end

%% floating base forward kinematics
[dim_fb, Xup, pos_idx] = fwd_kin_fb(model,q);
nb_pos = length(pos_idx);

%% Initialization
R_link  = cell(length(link_idx),1); 
switch class(q)
    case 'double'
        p_link  = zeros(nb_pos*length(link_idx),1);
    case 'casadi.SX'
        p_link  = casadi.SX.sym('pf',nb_pos*length(link_idx),1);
    case 'casadi.MX'
        p_link  = casadi.MX(zeros(nb_pos*length(link_idx),1));
    otherwise
        error('Invalid variable type for "q"')
end

X0  = cell(model.NB,1);
R0  = cell(model.NB,1);
p0  = cell(model.NB,1);




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
    indx = (nb_pos*i-(nb_pos-1)):(nb_pos*i);
    p_link(indx) = p0{link_idx(i)}(pos_idx);
    R_link{i} = R0{link_idx(i)};
end

X0_link = {X0{link_idx}};
