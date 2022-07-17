%% Calculate Derivatives of the Fourier Trajectory w.r.t. params

% Adapted from Bryan Dongik Lee, 2018
% Original code at: https://github.com/SNURobotics/optimal-excitation/tree/master/libraries/trajectory

% helper for testing arm regressor and regressor derivative functions

%% Inputs
% [Name]       [Description]                      [Size]
%  params       spline coefficients                2m*n    (parameter num * joints num)
%  w            base frequency (sin wt)            1*1
%  t            times to sample                    1*k     (1*desired number of frames)                 

%% Outputs
% [Name]       [Description]                                 [Size]
%  dq           spline pos derivative vector                  n*2mn*k
%  dqdot        (optional) spline vel derivative vector       n*2mn*k
%  dqddot       (optional) spline acc derivative vector       n*2mn*k       

%% Implementation
function [dq, dqdot, dqddot] = getFourierDerivative(params, w, t)
    n = size(params,2);
    m = floor(size(params,1)/2); % k = 1,2,...m
    
    num_t = size(t,2);
    
%     dq = zeros(n, 2*m, num_t);
%     dqdot = zeros(n, 2*m, num_t);
%     dqddot = zeros(n, 2*m, num_t);
    
    % for params as a 2mn long vector
    dq = zeros(n, 2*m*n, num_t);
    dqdot = zeros(n, 2*m*n, num_t);
    dqddot = zeros(n, 2*m*n, num_t);

    for i = 1:n
        for j = 1:num_t
            for k = 1:m
%                 dq(i,2*k-1,j)     =   sin(k*w*t(j));
%                 dq(i,2*k,j)       =   cos(k*w*t(j));
%                 dqdot(i,2*k-1,j)  =   k*w*cos(k*w*t(j));
%                 dqdot(i,2*k,j)    = - k*w*sin(k*w*t(j));
%                 dqddot(i,2*k-1,j) = - k*k*w*w*sin(k*w*t(j));                
%                 dqddot(i,2*k,j)   = - k*k*w*w*cos(k*w*t(j));  

                % for params as a 2mn long vector
                koff = (i-1)*2*m;
                dq(i,koff+2*k-1,j)     =   sin(k*w*t(j));
                dq(i,koff+2*k,j)       =   cos(k*w*t(j));
                dqdot(i,koff+2*k-1,j)  =   k*w*cos(k*w*t(j));
                dqdot(i,koff+2*k,j)    = - k*w*sin(k*w*t(j));
                dqddot(i,koff+2*k-1,j) = - k*k*w*w*sin(k*w*t(j));                
                dqddot(i,koff+2*k,j)   = - k*k*w*w*cos(k*w*t(j));  

            end
        end
    end
end