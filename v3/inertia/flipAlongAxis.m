function I = flipAlongAxis(I_in, axis)

%%% Get Pseudo Inertia
h = skew_spatial(I_in(1:3,4:6));
Ibar = I_in(1:3,1:3);
m = I_in(6,6);

if strcmp(class(I_in),'casadi.MX') % make sure this function works for optimization parameters
    P = casadi.MX.zeros(4,4);
    I = casadi.MX.eye(6);
else
    P = zeros(4,4);
    I = eye(6);
end

P(1:3,1:3) = 0.5*trace(Ibar)*eye(3) - Ibar;
P(1:3,4) = h;
P(4,1:3) = h';
P(4,4) = m;

%%% Flip Along Axis
X = eye(4);
if (axis == 'X')
    X(1, 1) = -1;
elseif (axis == 'Y')
    X(2, 2) = -1;
elseif (axis == 'Z')
    X(3, 3) = -1;
end
P = X * P * X;

m = P(4,4);
h = P(1:3,4);
E = P(1:3,1:3);
Ibar = trace(E) * eye(3) - E;
I(1:3,1:3) = Ibar;
I(1:3,4:6) = skew_spatial(h);
I(4:6,1:3) = skew_spatial(h)';
I(4:6,4:6) = m * eye(3);

end