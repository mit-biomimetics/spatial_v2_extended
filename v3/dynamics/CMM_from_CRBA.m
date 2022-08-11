function  [A] = CMM_from_CRBA( model, q)

% Algo from Wensing and Orin (IJHR)

% assert(strcmp(model.jtype(1),'Fb'), 'First joint should be floating base');


if ~isfield(model,'nq')
    model = postProcessModel(model);
end
qd = q(1)*0 + zeros(model.NV,1);
% if ~iscell(q)
%     [q,qd] = confVecToCell(model,q,qd);
% end



% Compute CMM from Mass Matrix

if strcmp(model.jtype(1),'Fb')
    % Joint model for floating base
    if ~iscell(q)
        [q,qd] = confVecToCell(model,q,qd);
    end
    
    [X10, Phi ] = jcalc( model.jtype{1}, q{1} );
    dim_fb = 6;
elseif strcmp(model.fb_type,'eul')
    % This is charles attempt to enable support for the way we make our models. It's almost right, but not quite. (see )
    [dim_fb, Xup, ~, R_world_to_body, rpy] = fwd_kin_fb(model,q);
    X10 = Xup{dim_fb};
 
    % B_euler_rate_to_omega_body = R_world_to_body*BmatF(rpy);
    B_euler_rate_to_omega_body = BmatF(rpy);
    B = R_world_to_body([1:3],:)*B_euler_rate_to_omega_body([1 2 3],:);
    Z33 = zeros(3,3);
    
    Phi = [Z33, B; R_world_to_body, Z33];


    if ~iscell(q)
        [q,qd] = confVecToCell(model,q,qd);
    end
    
else
    error('floating base type not supported')    
end

[H] = HandC(model,q,qd);

Psi = inv(Phi);
H11 = H(1:dim_fb,1:dim_fb);

IC  = Psi'* H11 * Psi;

[~, p1G] = mcI( IC ); % extract CoM position rel to FB

R1G = X10(1:3,1:3);
X1G = [R1G  zeros(3); skew(p1G)*R1G R1G];

A = X1G' * Psi' * H(1:dim_fb, :);