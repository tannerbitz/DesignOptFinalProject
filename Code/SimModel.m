classdef SimModel < handle
   properties
      % vehicle physical properties
      mass      % vehicle mass
      length    % vehicle length
      dm        % front/back load balance 
      Iz        % moment of inertia about z axis
      hcg       % center of gravity height
      wr        % wheel radius
      maxT      % maximum drivetrain torque
      mu        % road/tire friction coefficient
      C         % tire stiffness
      
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
       
       
       function update(obj,deltaT)
           % update vehicle states using ODE45
           
           States0 = [obj.X, obj.Y, obj.phi, obj.vx, obj.vy, obj.omegaB, obj.omegaW_fl, obj.omegaW_fr, obj.omegaW_rl, obj.omegaW_rr]';
           
           [~,newStates] = ode45(@(t,States) obj.calcRHS(t,States, delta,Torque), [0 deltaT], States0);
        
           [obj.X, obj.Y, obj.phi, obj.vx, obj.vy, obj.omegaB, obj.omegaW_fl, obj.omegaW_fr, obj.omegaW_rl, obj.omegaW_rr] = newStates(end,:) ;
       end
       
       function [States_dot] = calcRHS(obj,~,States, delta, torque)
           
           X_in = States(1);
           Y_in = States(2);
           phi_in = States(3);
           vx_in = States(4);
           vy_in = States(5);
           omegaB_in = States(6);
           omegaW_fl_in = States(7);
           omegaW_fr_in = States(8);
           omegaW_rl_in = States(9);
           omegaW_rr_in = States(10);
           
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
           [Ffl_long, Ffl_lat] = obj.getTireForces(omegaW_fl_in, vx_in, vy_in, Ffl_z, delta, torque);
           [Ffr_long, Ffr_lat] = obj.getTireForces(omegaW_fr_in, vx_in, vy_in, Ffr_z, delta, torque);
           [Frl_long, Frl_lat] = obj.getTireForces(omegaW_rl_in, vx_in, vy_in, Frl_z, 0, torque);
           [Frr_long, Frr_lat] = obj.getTireForces(omegaW_rr_in, vx_in, vy_in, Frr_z, 0, torque);
              
           % update Fx and Fy so it can be used in the next timestep;
           obj.Fx =  Frl_long + Frr_long - Ffl_lat*sin(delta) - Ffr_lat*sin(delta) ...
               + Ffl_long*cos(delta) + Ffr_long*cos(delta);
           
           obj.Fy = Frl_lat + Frr_lat + Ffl_lat*cos(delta) + Ffr_lat*cos(delta) ...
               + Ffl_long*sin(delta) + Ffr_long*sin(delta);
           
           % compute derivatives of vehicle states
           X_dot = vx_in*cos(phi_in) - vy_in*sin(phi_in);
           Y_dot = vx_in*sin(phi_in) + vy_in*cos(phi_in);
           phi_dot = omegaB_in;
           vx_dot  = 1/obj.mass*(obj.Fx - obj.mass*vy_in*obj.omegaB);
           vy_dot  = 1/obj.mass*(obj.Fy - obj.mass*vx_in*obj.omegaB);
           omegaB_dot = 1/obj.Iz*(Ffl_lat*lf*cos(delta) + Ffr_lat*lf*cos(delta) ...
               - Frl_lat*lr - Frr_lat*lr + (-Ffl_long - Frl_long + Ffr_long + Frr_long) * obj.width/2);
           omegaW_fl_dot = torque - obj.wr * Ffl_long;
           omegaW_fr_dot = torque - obj.wr * Ffr_long;
           omegaW_rl_dot = torque - obj.wr * Frl_long;
           omegaW_rr_dot = torque - obj.wr * Frr_long;
           
           States_dot = [X_dot;Y_dot; phi_dot; vx_dot; vy_dot; omegaB_dot; omegaW_fl_dot; omegaW_fr_dot; omegaW_rl_dot; omegaW_rr_dot];
           
       end
 
      function [F_long, F_lat] = getTireForces(obj, omegaW_in, vx_in, vy_in, Fz_in, delta_in, torque_in)
            %getTireForces calculates and returns the latitudinal and
            %longitudinal forces on a wheel.
            
            % Calculate sigmaLat/ sigmaLong/ sigma
            if sign(torque_in) > 0
                sigmaLong = (obj.wr*omegaW_in - vx_in)/(obj.wr*omegaW_in);
            else
                sigmaLong = (obj.wr*omegaW_in - vx_in)/(vx_in);
            end
            
            alpha = delta_in - atan(vy_in/vx_in);
            
            sigmaLat = vx_in/(obj.wr*omegaW_in)*tan(alpha);
            sigma = sqrt(sigmaLat^2 + sigmaLong^2);
            
            % Calculate theta
            theta = obj.C/(3*obj.mu*Fz_in);
            
            % Determine if tire is sliding/ total force
            thetaPoly = 3*theta*sigma - 3*(theta*sigma)^2 + (theta*sigma)^3;
            if thetaPoly < 1 %not sliding
                Ft = obj.mu*Fz_in*thetaPoly;
            else
                Ft = Fz;
            end
            
            % Calculate lat and long tire forces
            F_long = sigmaLong/sigma*Ft;
            F_lat = sigmaLat/sigma*Ft;
            
        end
       
       
   end
end