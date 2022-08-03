function [pf,Rf] = get_gc_position( model, q, gc)
% pf = get_gc_position( model, q, gc) computes the position "pf" (in the 
% world frame) of the ground contact point with index "gc" for the given
% robot "model" in the given joint configuration "q".
% The argument gc can be a list of ground contacts indexes. The pf will be
% a vector of size 3*number of ground contact indices given

%% Initialization of variables
X0  = cell(model.NB,1);
switch class(q)
    case 'double'
        pf  = zeros(3*length(gc),1);
        Rf  = zeros(3);
    case 'casadi.SX'
        pf  = casadi.SX.sym('pf',3*length(gc),1);
        Rf  = casadi.SX.sym('Rf',3,3);
    case 'casadi.MX'
        pf  = casadi.MX(zeros(3*length(gc),1));
        Rf  = casadi.MX(zeros(3,3));
    otherwise
        error('Invalid variable type for "q"')
end

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

for i = (dim_fb + 1):model.NB
    X0{i} = Xup{i} * X0{model.parent(i)}; % propagate xform from origin
end

for i = 1:length(gc)
    indx = (3*i-2):(3*i);
    [Rfi,pfi] = plux(model.gc_X{gc(i)} * X0{model.gc_parent(gc(i))}); % origin to foot translation, world coordinates
    pf(indx) = pfi;
end
Rf = Rfi;


