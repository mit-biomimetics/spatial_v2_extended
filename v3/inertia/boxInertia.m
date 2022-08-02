function I = boxInertia(mass, x)
    I = (norm(x)^2*eye(3) - diag(x.^2))*mass/12;
end