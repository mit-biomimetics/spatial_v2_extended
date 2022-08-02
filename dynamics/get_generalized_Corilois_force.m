function Cqd = get_generalized_Corilois_force(model, q, qd)
% Computes the generalized Coriolis forces (Cqd) in the inverse dynamics
% formulation (assuming the model treats its floating base as six 1DOF
% joints with euler angle convention for the orientation)
%
% @return Cqd (model.NB x 1) - generalized coriolis force (body frame)

if ~strcmp(model.fb_type,'eul')
    error('get_mass_matrix only works with euler angle-based floating base joint')
end

%% Forward Kinematics
R_world_to_body = rpyToRotMat(q(4:6))';
% for i = 1:5
%     Xup{i} = zeros(6,6);
%     v{i} = zeros(6,1);
%     c{i} = zeros(6,1);
% end
Xup{6} = [R_world_to_body zeros(3,3);...
    -R_world_to_body*skew_spatial(q(1:3)) R_world_to_body];
v{6} = qd(1:6); % [omega (body frame), vBody (body frame)

for i = 7:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
    vJ = S{i}*qd(i);
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    c{i} = crm(v{i}) * vJ;
end

%% Bias Accelerations
avp{6} = zeros(6,1);
for i = 7:model.NB
    avp{i} = Xup{i} * avp{model.parent(i)} + c{i};
end

%% Generalized Coriolis Force
switch class(q)
    case 'double'
        Cqd = zeros(model.NB,1);
    case 'casadi.SX'
        Cqd = casadi.SX.sym('Cqd',model.NB,1);
    case 'casadi.MX'
        Cqd = casadi.MX(zeros(model.NB,1));
    otherwise
        error('Invalid variable type for "q"')
end

% floating base force
Ifb = model.I{6};
hfb = Ifb*v{6};
fvp{6} = Ifb*avp{6} + crf(v{6})*hfb;

for i = 7:model.NB
   Ii = model.I{i};
   hi = Ii*v{i};
   fvp{i} = Ii*avp{i} + crf(v{i})*hi;
end

for i = model.NB:-1:7
   Cqd(i) = S{i}'*fvp{i};
   fvp{model.parent(i)} = fvp{model.parent(i)}  + Xup{i}'*fvp{i};
end

Cqd(1:6) = fvp{6};

