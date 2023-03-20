function  [H,Ic, IC] = get_mass_matrix( model, q)
% Computes mass matrix (H) in the inverse dynamics formulation (assuming
% the model treats its floating base as six 1DOF joints with euler angle 
% convention for the orientation)
%
% @return H (model.NB x model.NB) - mass matrix (body frame)
% @return Ic (6 x 6) - composite rigid body inertia (body frame)

if ~any(strcmp(model.fb_type,{'eul','planar'}))
    error('get_mass_matrix only works with euler angle-based floating base joint')
end

%% 2D or 3D
[Xup, S] = fwd_kin_fb(model, q);
dim_fb = model.NB_fb;

%% Forward Kinematics
for i = (dim_fb + 1):model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
end

%% Composite Inertia
IC = model.I;
for i = model.NB:-1:(dim_fb + 1)
    IC{model.parent(i)} = IC{model.parent(i)} + Xup{i}'*IC{i}*Xup{i};
end

%% Mass Matrix
switch class(q)
    case 'double'
        H = zeros(model.NB);
    case 'casadi.SX'
        H = casadi.SX.sym('H', model.NB, model.NB);
    case 'casadi.MX'
        H = casadi.MX(zeros(model.NB));
    otherwise
        error('Invalid variable type for "q"')
end


H(1:dim_fb,1:dim_fb) = S{dim_fb}' * IC{dim_fb} * S{dim_fb};

for i = (dim_fb + 1):model.NB
    fh = IC{i} * S{i};
    H(i,i) = S{i}' * fh;
    
    fh = Xup{i}' * fh;
    j = model.parent(i);
    while j > dim_fb
        H(i,j) = S{j}' * fh;
        H(j,i) = H(i,j);
        fh = Xup{j}' * fh;
        j = model.parent(j);
    end
    
    H(1:dim_fb,i) = S{dim_fb}' * fh;
    H(i,1:dim_fb) = fh' * S{dim_fb};
end

% Composite Rigid Body Inertia - Body frame
U = [eye(dim_fb) zeros(dim_fb, model.NB - dim_fb)];
Ic = U * H * U';
