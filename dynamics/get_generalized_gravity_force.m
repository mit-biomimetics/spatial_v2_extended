function G = get_generalized_gravity_force(model, q)
% Computes the generalized gravity forces (G) in the inverse dynamics
% formulation (assuming the model treats its floating base as six 1DOF
% joints with euler angle convention for the orientation)
%
% @return G (model.NB x 1) - generalized gravity force (body frame)

if ~strcmp(model.fb_type,'eul')
    error('get_mass_matrix only works with euler angle-based floating base joint')
end

%% Forward Kinematics
R_world_to_body = rpyToRotMat(q(4:6))';
% for i = 1:5
%     Xup{i} = zeros(6,6);
% end
Xup{6} = [R_world_to_body zeros(3,3);...
    -R_world_to_body*skew_spatial(q(1:3)) R_world_to_body];

for i = 7:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
end

%% Composite Inertia
IC = model.I;
for i = model.NB:-1:7
    IC{model.parent(i)} = IC{model.parent(i)} + Xup{i}'*IC{i}*Xup{i};
end

%% Generalized Gravity Force
switch class(q)
    case 'double'
        G = zeros(model.NB,1);
    case 'casadi.SX'
        G = casadi.SX.sym('G',model.NB,1);
    case 'casadi.MX'
        G = casadi.MX(zeros(model.NB,1));
    otherwise
        error('Invalid variable type for "q"')
end

ag{6} = Xup{6} * get_gravity(model);
G(1:6) = -IC{6} * ag{6};

for i = 7:model.NB
   ag{i} = Xup{i} * ag{model.parent(i)};
   G(i) = -S{i}'*(IC{i}*ag{i}); 
end

