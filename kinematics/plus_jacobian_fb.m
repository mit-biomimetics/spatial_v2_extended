function J = plus_jacobian_fb(q_fb)
    % This function returns the jacobian for addition on the non-Euclidean manifold that describes
    % the motion of the floating base
    %
    % q_plus = q + J * dq
    % 
    % where J is the plus jacobian and dq is an infinitesimal unit vector in the tangent space of 
    % the manifold containing q

    if (size(q_fb) == [6 1]) % 3D
        rpy = q_fb(4:6);
        R_body_to_world = rpyToRotMat(rpy);
        J = [zeros(3) R_body_to_world;
             Binv(rpy) * R_body_to_world zeros(3)];
    elseif (size(q_fb) == [3 1]) % Planar
        J = [[0; 0; 1] [rpyToRotMat(q_fb(3)); 0 0]];
    else
        error('Invalid dimension of floating base coordinates')
    end

end

function B = Binv(rpy)
    % Generates the inverse of B matrix where the B matrix is used to convert
    % euler rates psidot, thetadot and phidot into angular velocities in world
    % coordinates
    %
    % The inverse B matrix (calculated here) is used to convert from anuglar
    % velocity in world coordinates to euler rates
    %
    % Suffers from singularity at theta = +- pi/2
    % 
    % Equation: ThetaDot = B^(-1)*omega

    psi = rpy(3);
    theta = rpy(2);
    B = [cos(psi) / cos(theta) sin(psi) / cos(theta) 0;
         -sin(psi) cos(psi) 0;
         cos(psi) * tan(theta) sin(psi) * tan(theta) 1];
end
