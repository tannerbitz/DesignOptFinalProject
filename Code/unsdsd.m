F_lat = car.getTireForces(Massera)



function F_lat = getTireForcesMassera(obj, omegaW_in, delta_in, vy_in, vx_in, Fz_in, frontBackLength)
            %GETTIREFORCESMASSERA returns the lateral force on a given
            %wheel based on the Fiala Tire Model section of Massera (2015)
            % frontBackLength =  lr for front wheel
            % frontBackLength = -lf for back wheel
            alpha = atan2(vy_in+frontBackLength*omegaW_in, vx_in) - delta_in ;
            f_alpha = obj.C*tan(alpha);
            R_mu = obj.mu_k/obj.mu_s;
    
            if abs(alpha) <= atan(3*obj.mu_s*Fz_in/obj.C)
                F_lat = -f_alpha + (2-R_mu)/(3*obj.mu_s*Fz_in)*abs(f_alpha)*f_alpha ... 
                        - (1-2/3*R_mu)/(3*obj.mu_s*Fz_in)^2*f_alpha^3;
            else
                F_lat = -sign(alpha)*obj.mu_s*R_mu*Fz_in;
            end
        end