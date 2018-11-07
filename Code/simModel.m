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
      
      % vehicle states - U
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
      
   end
   
   methods
       
       function [U_dot] = getRHS(obj,control, Fx,Fy)
           T = control(1);      % commanded wheel torque
           delta = control(2);  % commanded steering angle  (radians)
           
           lf = obj.dm*obj.length;
           lr = (1-obj.dm)*obj.length;
           
           
           
           
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
           omegaW_fr_dot = T - obj.wr * Frl_long;
           omegaW_rr_dot = T - obj.wr * Frr_long;
           
           U_dot = [X_dot; Y_dot; phi_dot; vx_dot; vy_dot; omegaB_dot; omegaW_fl_dot; omegaW_fr_dot; omegaW_rr_dot];
           
       end
       
       
   end
end