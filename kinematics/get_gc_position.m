function [pf, Rf, X0, Xup] = get_gc_position( model, q, gc)
% pf = get_gc_position( model, q, gc) computes the position "pf" (in the 
% world frame) of the ground contact point with index "gc" for the given
% robot "model" in the given joint configuration "q".
% The argument gc can be a list of ground contacts indexes. The pf will be
% a vector of size 3*number of ground contact indices given

%% Initialization of variables
X0  = cell(model.NB,1);
Rf = cell(length(gc),1); 

if strcmp(model.fb_type,'eul')
    %% floating base forward kinematics
    [NB_fb, Xup, pos_idx] = fwd_kin_fb(model, q);
else
    NB_fb = 1;
    pos_idx = 1:3;
    [ XJ, ~ ] = jcalc( model.jtype{1}, q );
    Xup{1} = XJ * model.Xtree{1};
end

nb_pos = length(pos_idx);
switch class(q)
    case 'double'
        pf  = zeros(nb_pos*length(gc),1);
    case 'casadi.SX'
        pf  = casadi.SX.sym('pf',nb_pos*length(gc),1);
    case 'casadi.MX'
        pf  = casadi.MX(zeros(nb_pos*length(gc),1));
    otherwise
        error('Invalid variable type for "q"')
end


%% Joints Forward Kinematics
if ~iscell(q) || ~iscell(qd)
    [q] = confVecToCell(model,q);
end

for i = (NB_fb + 1):model.NB
    [ XJ, ~ ] = jcalc( model.jtype{i}, q{i} ); % HAVE TO USE [q_cell, v_cell, vd_cell, vd2_cell] = confVecToCell(model,q,v,vd, vd2)
    Xup{i} = XJ * model.Xtree{i};
end

X0{NB_fb} = Xup{NB_fb};

for i = (NB_fb + 1):model.NB
    X0{i} = Xup{i} * X0{model.parent(i)}; % propagate xform from origin
end

for i = 1:length(gc)
    indx = (nb_pos*i-(nb_pos-1)):(nb_pos*i);
    [Rfi,pfi] = plux(model.gc_X{gc(i)} * X0{model.gc_parent(gc(i))}); % origin to foot translation, world coordinates
    pf(indx) = pfi(pos_idx);
end
Rf = Rfi;
