classdef MPC < handle
   properties
      f
      A
      B
      Gamma_el
      Gamma_ec
      C_el
      C_ec
      
      States 
      Ts
      N     % Horizon length
       
   end
   
   methods
       function obj = MPC(~,Ts)
           addpath ./TaylorEqs/;
           
           obj.Ts = Ts;
       end
       
       function [] = setStates(obj,States, uN)
           N1 = obj.N;
           % set initial states to current sim states
           obj.States(:,1) = States;  
           
           % shift states by one sampling period
           for i = 2:N1-1
               obj.States(:,i) = obj.States(:,i+1);
           end
           
           % compute last state
           obj.States(:,N1) = obj.States(:,N1) + obj.Ts*obj.f(:,N1);
           
       end
       
       function [] = linearize(obj, car)
           
           N1 = obj.N;
           for i = 1:N1
               obj.Gamma_el =
               
               if alphaf < tan  & alphar <tan
               obj.fn(:,i) = 
               obj.A(:,:,i) =   f_fgrg(car.mass
               obj.B(:,:,i) =    
               
               
               elseif
                   
                   
               elseif
                   
               else
           end
       end
       
       function cost = cost(obj,inVec)
           
       end
       
       
       
       
   end
       
    
end
