function [pf, Rf, X0, Xup] = get_gc_position( model, q, gc)
% pf = get_gc_position( model, q, gc) computes the position "pf" (in the 
% world frame) of the ground contact point with index "gc" for the given
% robot "model" in the given joint configuration "q".
% The argument gc can be a list of ground contacts indexes. The pf will be
% a vector of size 3*number of ground contact indices given

%% Initialization of variables
X0  = cell(model.NB,1);
Rf = cell(length(gc),1); 

Xup = get_spatial_transforms(model, q);
nb_pos = length(model.pos_idx);

switch class(q)
    case 'double'
        pf  = zeros(nb_pos*length(gc),1);
    case 'casadi.SX'
        pf  = casadi.SX(nb_pos*length(gc),1);
    case 'casadi.MX'
        pf  = casadi.MX(nb_pos*length(gc), 1);
    otherwise
        error('Invalid variable type for "q"')
end

X0{model.NB_fb} = Xup{model.NB_fb};

for i = (model.NB_fb + 1):model.NB
    X0{i} = Xup{i} * X0{model.parent(i)}; % propagate xform from origin
end

for i = 1:length(gc)
    indx = (nb_pos*i-(nb_pos-1)):(nb_pos*i);
    [Rf{i}, pfi] = plux(model.gc_X{gc(i)} * X0{model.gc_parent(gc(i))}); % origin to foot translation, world coordinates
    pf(indx) = pfi(model.pos_idx);
end

if length(gc) == 1
    Rf = Rf{1};
end