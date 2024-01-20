classdef revolutePairDifferentialWithRotors
    %revoluteJointWithRotor 
    
    properties
        nv = 2
        nq = 2
        jointAxis
        rotorAxis
        
        gearRatio
        
        diffMat


        XtreeInternal
        output_body = 4
        bodies = 4
    end
    
    methods        
        function [Xup, S, Sd, v] = kinematics(obj, Xtree, q, qdot, vp)
            % Assuming that the generalized coordaintes are the realtive
            % angles. All the challenge comes from rotors.
            
            [XJ1 , S1 ] = jcalc( ['R' obj.jointAxis{1}] , q(1) );
            
            [XJ2 , S2 ] = jcalc( ['R' obj.jointAxis{2}] , q(2) );
           
            n1 = obj.gearRatio{1};
            n2 = obj.gearRatio{2};
            G = obj.diffMat; % 2x2 matrix that transforms differential output angles to differential input angles

            % get differential input angles from output angles
            % then multiply by respective gear ratios to get rotor angles
            qr1 = n1*G(1,:)*q; 
            qr2 = n2*G(2,:)*q;
            
            [XR1 , Sr1 ] = jcalc( ['R' obj.rotorAxis{1}], qr1 );
            [XR2 , Sr2 ] = jcalc( ['R' obj.rotorAxis{2}], qr2 );
            
            X1p  = XJ1 * Xtree(1:6,:);
            Xr1p = XR1 * Xtree(7:12,:);
            Xr2p = XR2 * Xtree(13:18,:);
            X21  = XJ2 * obj.XtreeInternal;
            
            Xup = [X1p ; Xr1p ; Xr2p ; X21*X1p];
            
            z = zeros(6,1);

            % coupling here is the hard part
            S  = [S1                          z; 
                  Sr1*n1*G(1,1)   Sr1*n1*G(1,2);
                  Sr2*n2*G(2,1)   Sr2*n2*G(2,2);
                  X21*S1                     S2
                  ];
              

            if nargout > 2
                % this is probably wrong, but I'm not using it
                v  = Xup*vp + S*qdot;
                v1 = v(1:6);
                vr1= v(7:12);
                vr2= v(13:18);
                v2 = v(19:24);
                
                S1d = crm(v1)*S1;
                S2d = crm(v2)*S2;
                Sr1d= crm(vr1)*Sr1;
                Sr2d= crm(vr2)*Sr2;
                
                Sd  = [ S1d                           z; 
                        Sr1d*n1*G(1,1)   Sr1d*n1*G(1,2);
                        Sr2d*n2*G(2,1)   Sr2d*n2*G(2,2);
                        X21*S1d                     S2d
                      ];              
            end
        end
        
    end
end

