function  Jdqd = BodyJdotqdot(model, q, qd, body_num, Xend)
    
    if ~isfield(model,'nq')
        model = postProcessModel(model);
    end

    if ~iscell(q) || ~iscell(qd)
        [q, qd] = confVecToCell(model,q, qd);
    end
    
    Jdqd = q{1}(1)*0 + zeros( size(Xend,1), 1 );

    for i = 1:model.NB
        [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
        vJ = S{i}*qd{i};
        Xup{i} = XJ * model.Xtree{i};
        if model.parent(i) == 0
          v{i} = vJ;
          avp{i} = zeros(6,1);
        else
          v{i} = Xup{i}*v{model.parent(i)} + vJ;
          avp{i} = Xup{i}*avp{model.parent(i)} + crm(v{i})*vJ;
        end
    end

    ac = Xend*avp{body_num};
    vc = Xend*v{body_num};
    Jdqd(1:3) = ac(1:3);
    Jdqd(4:6) = spatialToLinearAcceleration(ac, vc);
end

function acc = spatialToLinearAcceleration(a, v)
    acc = a(4:6) + cross(v(1:3),v(4:6));
end