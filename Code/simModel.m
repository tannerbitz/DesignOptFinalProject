classdef simModel < handle
   properties
      % vehicle physical properties
      mass      % vehicle mass
      length    % vehicle length
      dm        % front/back load balance 
      Iz        % moment of inertia about z axis
      hcg       % center of gravity height
      wr        % wheel radius
      maxT      % maximum wheel torque
      mu        % road/tire friction coefficient
      
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
       
       function [States_dot] = getRHS(obj,control)
           T = control(1);      % commanded wheel torque
           delta = control(2);  % commanded steering angle  (radians)
           
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
           [Ffl_long, Ffl_lat] = obj.getTireForces(obj.omegaW_fl,Ffl_z,T);
           [Ffr_long, Ffr_lat] = obj.getTireForces(obj.omegaW_fr,Ffr_z,T);
           [Frl_long, Frl_lat] = obj.getTireForces(obj.omegaW_rl,Frl_z,T);
           [Frr_long, Frr_lat] = obj.getTireForces(obj.omegaW_rr,Frr_z,T);
           
           % compute derivatives of vehicle states
           X_dot = obj.vx*cos(obj.phi) - obj.vy*sin(obj.phi);
           Y_dot = obj.vx*sin(obj.phi) + obj.vy*cos(obj.phi);
           phi_dot = obj.omegaB;
           vx_dot  = 1/obj.mass*(Frl_long + Frr_long - Ffl_lat*sin(delta) - Ffr_lat*sin(delta) ...
               + Ffl_long*cos(delta) + Ffr_long*cos(delta) - obj.mass*obj.vy*obj.omegaB);
           vy_dot  = 1/obj.mass*(Frl_lat + Frr_lat + Ffl_lat*sin(delta) + Ffr_lat*sin(delta) ...
               + Ffl_long*cos(delta) + Ffr_long*cos(delta) - obj.mass*obj.vy*obj.omegaB);
           omegaB_dot = 1/obj.Iz*(Ffl_lat*lf*cos(delta) + Ffr_lat*lf*cos(delta) ...
               - Frl_lat*lr - Frr_lat*lr + (-Ffl_long - Frl_long ...
               + Ffr_long + Frr_long) * obj.width/2);
           omegaW_fl_dot = T - obj.wr * Ffl_long;
           omegaW_fr_dot = T - obj.wr * Ffr_long;
           omegaW_rl_dot = T - obj.wr * Frl_long;
           omegaW_rr_dot = T - obj.wr * Frr_long;
           
           States_dot = [X_dot; Y_dot; phi_dot; vx_dot; vy_dot; omegaB_dot; omegaW_fl_dot; omegaW_fr_dot; omegaW_rl_dot; omegaW_rr_dot];
           
       end
       
      function  [F_long, F_lat] = getTireForces(obj, omegaW, Fz, delta, T)
            %getTireForces calculates and returns the latitudinal and
            %longitudinal forces on a wheel. The inputs are the normal
            %force on the tire, Fz, and a flag that specifies whether the
            %car is accelerating or braking. 
            
            % Calculate sigmaLat/ sigmaLong/ sigma
            if sign(T) > 0
                sigmaLong = (obj.wr*omegaW - obj.vx)/(obj.wr*obj.omegaW);
            else
                sigmaLong = (obj.wr*omegaW - obj.vx)/(obj.vx);
            end
            
            alpha = delta - atan(obj.vy/obj.vx);
            
            sigmaLat = obj.vx/(obj.wr*omegaW)*tan(alpha);
            sigma = sqrt(sigmaLat^2 + sigmaLong^2);
            
            % Calculate theta
            theta = obj.C/(3*obj.mu*Fz);
            
            % Determine if tire is sliding/ total force
            thetaPoly = 3*theta*sigma - 3*(theta*sigma)^2 + (theta*sigma)^3;
            if thetaPoly < 1 %not sliding
                Ft = obj.mu*Fz*thetaPoly;
            else
                Ft = Fz;
            end
            
            % Calculate lat and long tire forces
            F_long = sigmaLong/sigma*Ft;
            F_lat = sigmaLat/sigma*Ft;
            
        end
       
       
   end
end