function I = spatialInertia(a, b, c)
% Returns the spatial inertia of a rigid body.
% In the case of 1 argument, a is the 10 inertia parameters that
% uniquely describe the links inertia properties
% We currently ignore off diagonal terms in the body inertia

% In the case of 3 arguments, a = mass, b = com, c = 3x3 body inertia

if nargin == 1
    I = eye(6);
    I(1,1) = a(5);
    I(2,2) = a(6);
    I(3,3) = a(7);
    % I(1:3,1:3) = [a(5) a(10) a(9);
    %     a(10) a(6) a(8);
    %     a(9) a(8) a(7)];
    
    cSkew = skew_spatial([a(2),a(3),a(4)]);
    I(1:3,4:6) = cSkew;
    I(4:6,1:3) = -cSkew';
    
    I(4:6,4:6) = a(1)*eye(3);
    
else
    
    cSkew = skew_spatial(b);
    I = [c + a*(cSkew*cSkew'), a*cSkew;
        a*cSkew', a * eye(3)];
    
end