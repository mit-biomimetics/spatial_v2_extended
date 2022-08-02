function [pf,Rf] = get_body_relative_gc_position( model, qj, gc)
% Computes the position "pf" of the ground contact point with index "gc"
% for the given robot "model" in the given joint configuration "q".
%
% @param model, spatial v2 floating base model
% @param qj, vector of joint positions for the robot (no floating base)
% @param gc, index of the ground contact of interest
% 
% @return pf, body relative position of the gc-th ground contact point
% (expressed in the body frame)
% @return Rf, rotation matrix of the ground contact frame relative to the
% body frame

if ~strcmp(model.fb_type,'eul')
    error('get_gc_position only works with euler angle-based floating base joint')
end

%%
Xup = cell(model.NB,1);
X0  = cell(model.NB,1);
switch class(qj)
    case 'double'
        pf  = zeros(3,1);
        Rf  = zeros(3);
    case 'casadi.SX'
        pf  = casadi.SX.sym('pf',3,1);
        Rf  = casadi.SX.sym('Rf',3,3);
    case 'casadi.MX'
        pf  = casadi.MX(zeros(3,1));
        Rf  = casadi.MX(zeros(3,3));
    otherwise
        error('Invalid variable type for "qj"')
end

%% Forward Kinematics
% for i = 1:5
%     Xup{i} = zeros(6,6);
% end
Xup{6} = [eye(3) zeros(3,3);...
    zeros(3,3) eye(3)];

for i = 7:model.NB
    [ XJ, ~ ] = jcalc( model.jtype{i}, qj(i-6) );
    Xup{i} = XJ * model.Xtree{i};
end

X0{6} = Xup{6};
for i = 7:model.NB
    X0{i} = Xup{i} * X0{model.parent(i)}; % propagate xform from origin
end

[Rf,pf] = plux(model.gc_X{gc} * X0{model.gc_parent(gc)});
