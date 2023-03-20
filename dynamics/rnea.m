function tau = rnea(model, q, qd, qdd, f_ext)

    % ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
    % ID(model,q,qd,qdd,f_ext) calculates the inverse dynamics of a kinematic
    % tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
    % of joint position, velocity and acceleration variables; and the return
    % value is a vector of joint force variables.  f_ext is an optional
    % argument specifying the external forces acting on the bodies.  It can be
    % omitted if there are no external forces.  The format of f_ext is
    % explained in the source code of apply_external_forces.

    if ~any(strcmp(model.fb_type,{'eul','planar'}))
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

    %% 2D or 3D
    [Xup, S] = fwd_kin_fb(model, q);
    Xa = Xup;
    dim_fb = model.NB_fb;

    %% Forward Kinematics
    a_grav = [0 0 0 model.gravity]';

    v{dim_fb} = S{dim_fb} * qd(1:dim_fb);
    a{dim_fb} = Xup{dim_fb} * (-a_grav) + S{dim_fb} * qdd(1:dim_fb);
    h{dim_fb} = model.I{dim_fb} * v{dim_fb};
    f{dim_fb} = model.I{dim_fb} * a{dim_fb} + crf(v{dim_fb}) * h{dim_fb};

    for i = (dim_fb + 1):model.NB
        [XJ, S{i}] = jcalc(model.jtype{i}, q(i));
        vJ = S{i} * qd(i);
        Xup{i} = XJ * model.Xtree{i};
        Xa{i} = Xup{i} * Xa{model.parent(i)};
        v{i} = Xup{i} * v{model.parent(i)} + vJ;
        a{i} = Xup{i} * a{model.parent(i)} + S{i} * qdd(i) + crm(v{i}) * vJ;
        h{i} = model.I{i} * v{i};
        f{i} = model.I{i} * a{i} + crf(v{i}) * h{i};
    end

    %% Apply External Forces
    if nargin == 5
        if length(f_ext) > 0
            for i = dim_fb:model.NB

                if length(f_ext{i}) > 0
                    Xa_force = [Xa{i}(1:3, 1:3) Xa{i}(4:6, 1:3); zeros(3, 3) Xa{i}(4:6, 4:6)];
                    f{i} = f{i} - Xa_force * f_ext{i};
                end
            end
        end
    end

    %% Backward Pass
    for i = model.NB:-1:(dim_fb + 1)
        p = model.parent(i);
        tau(i) = S{i}.' * f{i};
        f{p} = f{p} + Xup{i}.' * f{i};
    end

    tau(1:dim_fb) = S{dim_fb}.' * f{dim_fb};

end
