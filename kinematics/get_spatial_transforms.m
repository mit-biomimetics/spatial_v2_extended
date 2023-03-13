function Xup = get_spatial_transforms(model, q)

    %% Initialization
    if ~strcmp(model.fb_type, 'eul')
        error('get_gc_position only works with euler angle-based floating base joint')
    end

    %% floating base forward kinematics
    [dim_fb, Xup, pos_idx] = fwd_kin_fb(model, q);
    nb_pos = length(pos_idx);

    %% Joints Forward Kinematics
    for i = (dim_fb + 1):model.NB
        [XJ, ~] = jcalc(model.jtype{i}, q(i));
        Xup{i} = XJ * model.Xtree{i};
    end

end
