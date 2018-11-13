classdef SimModel < handle
    properties
        % vehicle physical properties
        mass      % vehicle mass
        length    % vehicle length (wheelbase)
        width     % car width ( between wheels)
        dm        % front/back load balance
        Iz        % moment of inertia about z axis
        hcg       % center of gravity height
        wr        % wheel radius
        Iw        % tire/wheel/drive train rotational intertia
        mu_k      % kinetic road/tire friction coefficient
        mu_s      % static road/tire friction coefficient
        C         % tire stiffness
        F_long_accel   % max (positive) acceleration force on each wheel (N)
        F_long_brake   % max (negative) braking force on each wheel (N)
        delta_max      % max steering angle
        delta_min      % min steering angle
        v_max          % max velocity
        v_min          % min velocity
        refA      % aerodynamics reference area
        Cd        % drag coefficient
        rhoAir    % density of air
        
        % vehicle states - States
        X         % global X coordinate
        Y         % global Y coordinate
        phi       % vehicle heading relative to global X axis  [rad]
        vx        % body x velocity
        vy        % body y velocity
        omegaB    % body rotation rate
        omegaW_fl % front left wheel rotation rate
        omegaW_fr % front right wheel rotation rate
        omegaW_rl % front left wheel rotation rate
        omegaW_rr % front right wheel rotation rate
        
        Fx        % total x-force on vehicle from previous interation
        Fy        % total y-force on vehicle from previous interation
        
    end
    
    methods
        
        function obj = SimModel(obj)
            
        end
        
        
        function update(obj,deltaT, delta, F_long)
            % update vehicle states using ODE45
            
            States0 = [obj.X, obj.Y, obj.phi, obj.vx, obj.vy, obj.omegaB]';
            
            
            [~,newStates] = ode45(@(t,States) obj.calcRHS(t,States, delta,F_long), [0 deltaT], States0);
            
            obj.X = newStates(end,1);
            obj.Y = newStates(end,2);
            obj.phi = newStates(end,3);
            obj.vx = newStates(end,4);
            obj.vy = newStates(end,5);
            obj.omegaB = newStates(end,6);
        end
        
        function [States_dot] = calcRHS(obj,~,States, delta, F_long)
            
            X_in = States(1);
            Y_in = States(2);
            phi_in = States(3);
            vx_in = States(4);
            vy_in = States(5);
            omegaB_in = States(6);

            
            % compute lengths from cg to front and rear end
            lf = obj.dm*obj.length;
            lr = (1-obj.dm)*obj.length;
            
            % compute tire normal force using totoal forces from previous
            % timestep
            Ffl_z = .5*obj.mass*9.81*obj.dm - obj.Fx*obj.hcg/obj.length - obj.Fy*obj.hcg/obj.width;
            Ffr_z = .5*obj.mass*9.81*obj.dm - obj.Fx*obj.hcg/obj.length + obj.Fy*obj.hcg/obj.width;
            Frl_z = .5*obj.mass*9.81*(1-obj.dm) + obj.Fx*obj.hcg/obj.length - obj.Fy*obj.hcg/obj.width;
            Frr_z = .5*obj.mass*9.81*(1-obj.dm) + obj.Fx*obj.hcg/obj.length + obj.Fy*obj.hcg/obj.width;
            
            % compute tire longitudenal and lateral forces with Fiala tire
            % model
            Ffl_lat = obj.getTireForcesMassera(omegaB_in, delta, vy_in, vx_in, Ffl_z, lr);
            Ffr_lat = obj.getTireForcesMassera(omegaB_in, delta, vy_in, vx_in, Ffr_z, lr);
            Frl_lat = obj.getTireForcesMassera(omegaB_in, 0, vy_in, vx_in, Frl_z, -lf);
            Frr_lat = obj.getTireForcesMassera(omegaB_in, 0, vy_in, vx_in, Frr_z, -lf);

            Ffl_long = F_long;
            Ffr_long = F_long;
            Frl_long = F_long;
            Frr_long = F_long;


            % update Fx and Fy so it can be used in the next timestep;
            obj.Fx =  Frl_long + Frr_long - Ffl_lat*sin(delta) - Ffr_lat*sin(delta) ...
                + Ffl_long*cos(delta) + Ffr_long*cos(delta);
            
            obj.Fy = Frl_lat + Frr_lat + Ffl_lat*cos(delta) + Ffr_lat*cos(delta) ...
                + Ffl_long*sin(delta) + Ffr_long*sin(delta);
            
            % compute derivatives of vehicle states
            X_dot = vx_in*cos(phi_in) - vy_in*sin(phi_in);
            Y_dot = vx_in*sin(phi_in) + vy_in*cos(phi_in);
            phi_dot = omegaB_in;
            vx_dot  = 1/obj.mass*(obj.Fx + obj.mass*vy_in*obj.omegaB - sign(vx_in)*.5*obj.rhoAir*obj.Cd*obj.refA*vx_in^2);
            vy_dot  = 1/obj.mass*(obj.Fy - obj.mass*vx_in*obj.omegaB);
            omegaB_dot = 1/obj.Iz*(Ffl_lat*lf*cos(delta) + Ffr_lat*lf*cos(delta) ...
                - Frl_lat*lr - Frr_lat*lr + (-Ffl_long*sin(delta) - Frl_long + Ffr_long*sin(delta) + Frr_long) * obj.width/2);

            
            States_dot = [X_dot,Y_dot, phi_dot, vx_dot, vy_dot, omegaB_dot]';
            
        end
        
        function [F_long, F_lat] = getTireForces(obj, omegaW_in, vx_in, vy_in, Fz_in, delta_in, torque_in)
            %getTireForces calculates and returns the latitudinal and
            %longitudinal forces on a wheel.
            eps = 10^-13;
            % Calculate sigmaLat/ sigmaLong/ sigma
            if sign(torque_in) >= 0
                sigmaLong = (obj.wr*omegaW_in - vx_in)/(obj.wr*omegaW_in+eps);
            else
                sigmaLong = (obj.wr*omegaW_in - vx_in)/(vx_in+eps);
            end
            
            
            alpha = delta_in - atan2(vy_in,vx_in);
            
            sigmaLat = vx_in/(obj.wr*omegaW_in+eps)*tan(alpha);
            sigma = sqrt(sigmaLat^2 + sigmaLong^2);
            
            % Calculate theta
            theta = obj.C/(3*obj.mu*Fz_in);
            
            % Determine if tire is sliding/ total force
            thetaPoly = 3*theta*sigma - 3*(theta*sigma)^2 + (theta*sigma)^3;
            if sigma < .15 %thetaPoly < 1 %not sliding
                Ft = obj.mu*Fz_in*thetaPoly;
            else
                Ft = obj.mu*Fz_in;
            end
            
            % Calculate lat and long tire forces
            F_long = sigmaLong/(sigma+eps)*Ft;
            F_lat = sigmaLat/(sigma+eps)*Ft;
            
        end
        

        function F_lat = getTireForcesMassera(obj, omegaW_in, delta_in, vy_in, vx_in, Fz_in, frontBackLength)
            %GETTIREFORCESMASSERA returns the lateral force on a given
            %wheel based on the Fiala Tire Model section of Massera (2015)
            % frontBackLength =  lr for front wheel
            % frontBackLength = -lf for back wheel
            alpha = atan2(vy_in+frontBackLength*omegaW_in, vx_in) - delta_in ;
            f_alpha = obj.C*tan(alpha);
            R_mu = obj.mu_k/obj.mu_s;
    
            if abs(alpha) <= atan(3*obj.mu_s*Fz_in/obj.C)
                F_lat = -f_alpha + (2-R_mu)/(3*obj.mu_s*Fz_in)*abs(f_alpha)*f_alpha - ...
                        (1-2/3*R_mu)/(3*obj.mu_s*Fz_in)^2*f_alpha^3;
            else
                F_lat = -sign(alpha)*obj.mu_s*Fz_in;
            end



        end


        function [] = plotCar(obj)
            deltax1 = [obj.length/2, -obj.length/2, -obj.length/2, obj.length/2, obj.length/2+obj.length/4];
            deltay1 = [obj.width/2, obj.width/2, -obj.width/2, -obj.width/2, 0];
            
            deltax = cos(obj.phi)*deltax1 - sin(obj.phi)*deltay1;
            deltay = sin(obj.phi)*deltax1 + cos(obj.phi)*deltay1;
            
            x = obj.X + deltax;
            y = obj.Y + deltay;
            
            patch(x,y,sqrt(obj.vx^2 + obj.vy^2))
            colormap(jet)
        end
        
        
    end
end
