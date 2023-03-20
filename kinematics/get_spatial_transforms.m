function Xup = get_spatial_transforms(model, q)

    %% Initialization

    if any(strcmp(model.fb_type,{'eul','planar'}))
        %% floating base forward kinematics
        [Xup] = fwd_kin_fb(model, q);
    else
        [ XJ, ~ ] = jcalc( model.jtype{1}, q );
        Xup{1} = XJ * model.Xtree{1};
    end   
    
    %% Joints Forward Kinematics
    if ~iscell(q) || ~iscell(qd)
        [q] = confVecToCell(model,q);
    end
    
    for i = (model.NB_fb + 1):model.NB
        [ XJ, ~ ] = jcalc( model.jtype{i}, q{i} );
        Xup{i} = XJ * model.Xtree{i};
    end
end
