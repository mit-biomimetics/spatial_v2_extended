function  [H,Ic, IC] = get_mass_matrix( model, q)
% Computes mass matrix (H) in the inverse dynamics formulation (assuming
% the model treats its floating base as six 1DOF joints with euler angle 
% convention for the orientation)
%
% @return H (model.NB x model.NB) - mass matrix (body frame)
% @return Ic (6 x 6) - composite rigid body inertia (body frame)

if ~strcmp(model.fb_type,'eul')
    error('get_mass_matrix only works with euler angle-based floating base joint')
end

%% Forward Kinematics
R_world_to_body = rpyToRotMat(q(4:6))';
for i = 1:5
    Xup{i} = zeros(6,6);
end

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

%% Mass Matrix
switch class(q)
    case 'double'
        H = zeros(model.NB);
    case 'casadi.SX'
        H = casadi.SX.sym('H',model.NB,model.NB);
    case 'casadi.MX'
        H = casadi.MX(zeros(model.NB));
    otherwise
        error('Invalid variable type for "q"')
end

H(1:6,1:6) = IC{6};

for i = 7:model.NB
    fh = IC{i} * S{i};
    H(i,i) = S{i}' * fh;
    
    fh = Xup{i}' * fh;
    j = model.parent(i);
    while j > 6
        H(i,j) = S{j}' * fh;
        H(j,i) = H(i,j);
        fh = Xup{j}' * fh;
        j = model.parent(j);
    end
    
    H(1:6,i) = fh;
    H(i,1:6) = fh';
end

% Composite Rigid Body Inertia - Body frame
U = [eye(6) zeros(6,model.NB-6)];
Ic = U*H*U';
