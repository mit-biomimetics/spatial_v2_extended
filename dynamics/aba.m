function qdd = aba(model, q, qd, tau, f_ext)

    if ~any(strcmp(model.fb_type,{'eul','planar'}))
        error('aba only works with euler angle-based floating base joint')
    end

    switch class(q)
        case 'double'
            qdd = zeros(model.NB, 1);
        case 'casadi.SX'
            qdd = casadi.SX.sym('qdd', model.NB, 1);
        case 'casadi.MX'
            qdd = casadi.MX(zeros(model.NB, 1));
        otherwise
            error('Invalid variable type for "q"')
    end

    symzero = q(1) * 0;

    %% 2D or 3D
    [Xup, S] = fwd_kin_fb(model, q);
    Xa = Xup;
    dim_fb = model.NB_fb;

    %% Forward Kinematics
    a_grav = [0 0 0 model.gravity]';

    v{dim_fb} = S{dim_fb} * qd(1:dim_fb);
    c{dim_fb} = symzero * v{dim_fb};
    IA{dim_fb} = symzero + model.I{dim_fb};
    pA{dim_fb} = crf(v{dim_fb}) * model.I{dim_fb} * v{dim_fb};

    for i = (dim_fb + 1):model.NB
        [XJ, S{i}] = jcalc(model.jtype{i}, q(i));
        vJ = S{i} * qd(i);
        Xup{i} = XJ * model.Xtree{i};
        Xa{i} = Xup{i} * Xa{model.parent(i)};
        v{i} = Xup{i} * v{model.parent(i)} + vJ;
        c{i} = crm(v{i}) * vJ;
        IA{i} = symzero + model.I{i};
        pA{i} = crf(v{i}) * model.I{i} * v{i};
    end

    %% Apply External Forces
    if nargin == 5
        if length(f_ext) > 0
            for i = dim_fb:model.NB

                if length(f_ext{i}) > 0
                    Xa_force = [Xa{i}(1:3, 1:3) Xa{i}(4:6, 1:3); zeros(3, 3) Xa{i}(4:6, 4:6)];
                    pA{i} = pA{i} - Xa_force * f_ext{i};
                end
            end
        end
    end

    %% Backward Pass
    for i = model.NB:-1:(dim_fb + 1)
        p = model.parent(i);

        U{i} = IA{i} * S{i};
        d{i} = S{i}.' * U{i};
        u{i} = tau(i) - S{i}.' * pA{i} - U{i}.' * c{i};
        U{i} = Xup{i}.' * U{i};

        Ia = Xup{i}.' * IA{i} * Xup{i};
        Ia = Ia - U{i} / d{i} * U{i}.';

        pa = Xup{i}.' * (pA{i} + IA{i} * c{i});
        pa = pa + U{i} * (d{i} \ u{i});

        IA{p} = IA{p} + Ia;
        pA{p} = pA{p} + pa;

    end

    U{dim_fb} = IA{dim_fb} * S{dim_fb};
    d{dim_fb} = S{dim_fb}.' * U{dim_fb};
    u{dim_fb} = tau(1:dim_fb) - S{dim_fb}.' * pA{dim_fb} - U{dim_fb}.' * c{dim_fb};
    U{dim_fb} = Xup{dim_fb}.' * U{dim_fb};

    %% Forward Pass
    ap = -a_grav;
    qdd(1:dim_fb) = d{dim_fb} \ (u{dim_fb} - U{dim_fb}.' * ap);
    a{dim_fb} = Xup{dim_fb} * ap + S{dim_fb} * qdd(1:dim_fb) + c{dim_fb};

    for i = (dim_fb + 1):model.NB
        p = model.parent(i);
        ap = a{p};
        qdd(i) = d{i} \ (u{i} - U{i}.' * ap);
        a{i} = Xup{i} * ap + S{i} * qdd(i) + c{i};
    end

end
