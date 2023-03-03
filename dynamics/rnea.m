function tau = rnea(model, q, qd, qdd, f_ext)

    % ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
    % ID(model,q,qd,qdd,f_ext) calculates the inverse dynamics of a kinematic
    % tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
    % of joint position, velocity and acceleration variables; and the return
    % value is a vector of joint force variables.  f_ext is an optional
    % argument specifying the external forces acting on the bodies.  It can be
    % omitted if there are no external forces.  The format of f_ext is
    % explained in the source code of apply_external_forces.

    if ~strcmp(model.fb_type, 'eul')
        error('get_mass_matrix only works with euler angle-based floating base joint')
    end

    switch class(q)
        case 'double'
            tau = zeros(model.NB, 1);
        case 'casadi.SX'
            tau = casadi.SX.sym('tau', model.NB, 1);
        case 'casadi.MX'
            tau = casadi.MX(zeros(model.NB, 1));
        otherwise
            error('Invalid variable type for "q"')
    end

    %% Forward Kinematics
    a_grav = [0 0 0 model.gravity]';
    R_world_to_body = rpyToRotMat(q(4:6))';

    for i = 1:5
        Xup{i} = zeros(6, 6);
        v{i} = zeros(6, 1);
        a{i} = zeros(6, 1);
        h{i} = zeros(6, 1);
        f{i} = zeros(6, 1);
    end

    S{6} = eye(6);
    Xup{6} = [R_world_to_body zeros(3, 3); ...
                  -R_world_to_body * skew_spatial(q(1:3)) R_world_to_body];
    v{6} = S{6} * qd(1:6);
    a{6} = Xup{6} * (-a_grav) + S{6} * qdd(1:6);
    h{6} = model.I{6} * v{6};
    f{6} = model.I{6} * a{6} + crf(v{6}) * h{6};

    for i = 7:model.NB
        [XJ, S{i}] = jcalc(model.jtype{i}, q(i));
        vJ = S{i} * qd(i);
        Xup{i} = XJ * model.Xtree{i};
        v{i} = Xup{i} * v{model.parent(i)} + vJ;
        a{i} = Xup{i} * a{model.parent(i)} + S{i} * qdd(i) + crm(v{i}) * vJ;
        h{i} = model.I{i} * v{i};
        f{i} = model.I{i} * a{i} + crf(v{i}) * h{i};
    end

    %% Backward Pass
    for i = model.NB:-1:7
        p = model.parent(i);
        tau(i) = S{i}.' * f{i};
        f{p} = f{p} + Xup{i}.' * f{i};
    end

    tau(1:6) = S{6}.' * f{6};

end
