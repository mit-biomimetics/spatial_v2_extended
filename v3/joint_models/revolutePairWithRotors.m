classdef revolutePairWithRotors
    %revoluteJointWithRotor 
    
    properties
        nv = 2
        nq = 2
        jointAxis
        rotorAxis
        
        gearRatio
        beltRatio

        XtreeInternal
        output_body = 4
        bodies = 4
    end
    
    methods        
        function [Xup, S, Sd, v] = kinematics(obj, Xtree, q, qdot, vp)
            % Assuming that the generalized coordinates are the relative
            % angles. All the challenge comes from rotors.
            
            [XJ1 , S1 ] = jcalc( ['R' obj.jointAxis{1}] , q(1) );
            
            [XJ2 , S2 ] = jcalc( ['R' obj.jointAxis{2}] , q(2) );
           
            gr1 = obj.gearRatio{1};
            gr2 = obj.gearRatio{2};
            b1_1  = obj.beltRatio{1};
            b2_1  = obj.beltRatio{2};
            b2_2  = obj.beltRatio{3};
            
            G = zeros(2,2);
            G(1,1) = gr1*b1_1;
            G(2,1) = gr2*b2_1;
            G(2,2) = gr2*b2_1*b2_2;
            
            qr1 = G(1,:)*q;
            qr2 = G(2,:)*q;
            
            [XR1 , Sr1 ] = jcalc( ['R' obj.rotorAxis{1}], qr1 );
            [XR2 , Sr2 ] = jcalc( ['R' obj.rotorAxis{2}], qr2 );
            
            X1p  = XJ1 * Xtree(1:6,:);
            Xr1p = XR1 * Xtree(7:12,:);
            Xr2p = XR2 * Xtree(13:18,:);
            X21  = XJ2 * obj.XtreeInternal;
            
            Xup = [X1p ; Xr1p ; Xr2p ; X21*X1p];
            
            z = zeros(6,1);
            S  = [S1            z; 
                  Sr1*G(1,1)    z;
                  Sr2*G(2,1)    Sr2*G(2,2);
                  X21*S1        S2
                  ];
              
            if nargout > 2
                v  = Xup*vp + S*qdot;
                v1 = v(1:6);
                vr1= v(7:12);
                vr2= v(13:18);
                v2 = v(19:24);
                
                S1d = crm(v1)*S1;
                S2d = crm(v2)*S2;
                Sr1d= crm(vr1)*Sr1;
                Sr2d= crm(vr2)*Sr2;
                
                Sd  = [ S1d           z; 
                        Sr1d*G(1,1)   z;
                        Sr2d*G(2,1)   Sr2d*G(2,2);
                        X21*S1d       S2d
                      ];              
            end
        end
        
    end
end

