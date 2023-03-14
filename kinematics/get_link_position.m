function [p_link,R_link,X0_link] = get_link_position( model, q, link_idx)
% p_link = get_link_position( model, q, gc) computes the position "p_link"
% (in the world frame) of the origin of frame corresponding to the 
% robot's link with index "link_idx" for the given robot "model" in 
% the given joint configuration "q".

Xup = get_spatial_transforms(model, q);


X0{model.fb_dim} = Xup{model.fb_dim};

for i = (model.fb_dim + 1):model.NB
    X0{i} = Xup{i} * X0{model.parent(i)};                      % propagate xform from origin
    [R0{i}, p0{i}] = plux(X0{i});                                % rotation from origin, translation from origin
end

nb_pos = length(model.pos_idx);
for i = 1:length(link_idx)
    indx = (nb_pos*i-(nb_pos-1)):(nb_pos*i);
    p_link(indx) = p0{link_idx(i)}(model.pos_idx);
    R_link{i} = R0{link_idx(i)};
end

X0_link = {X0{link_idx}};
