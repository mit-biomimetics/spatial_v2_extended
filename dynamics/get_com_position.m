function com = get_com_position(model, q);

[~, Ic] = get_mass_matrix(model, q);

if (size(Ic) == [3, 3])
    R_body_to_world = [cos(q(3)) -sin(q(3)); sin(q(3)) cos(q(3))]';
    com_offset = [Ic(1, 2); Ic(2, 3)] ./ model.mass;
    com = q(1:2) + R_body_to_world * com_offset;
elseif (size(Ic) == [6, 6])
    R_body_to_world = rpyToRotMat(q(4:6));
    com_offset = [Ic(3, 5); Ic(1, 6); Ic(2, 4)] ./ model.mass;
    com = q(1:3) + R_body_to_world * com_offset;
else
    error('Invalid size for composite inertia of floating base')
end

