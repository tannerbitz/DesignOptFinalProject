classdef MPC < handle
   properties
      f
      A
      B
      Gamma_el
      Gamma_ec
      C_el
      C_ec
      
      States    % (6xN) [X, Y, varphi, vx, vy, omegaB]'
      thetaA    % (1xN) [thetaA]
      U         % (2xN) [Flong, delta]'
      
      Ts
      N         % Horizon length
       
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
            
            Fz_f = car.mass*9.81*car.dm;
            Fz_r = car.mass*9.81*(1-car.dm);

            N1 = obj.N;
            for i = 1:N1
                alphaf = atan2(obj.States(5, i) + car.lr*obj.States(6, i), obj.States(4, i)) - obj.U(2, i);
                alphar = atan2(obj.States(5, i) - car.lf*obj.States(6, i), obj.States(4, i));
    
                R = car.mu_k/car.mu_s;
                alphaFMax = atan2(3*car.mu_s*Fz_f, car.C);
                alphaRMax = atan2(3*car.mu_s*Fz_r, car.C);

                %%%   INPUT VECTOR %%%%%
                %    C,Flong,Fz_f,Fz_r,Iz,R,delta,lf,lr,m,mu,omegaB,varphi,vx,vy
                C_in = car.C;
                Flong_in = obj.U(1, i);
                Fz_f_in = Fz_f;
                Fz_r_in = Fz_r;
                Iz_in = car.Iz;
                R_in = R;
                delta_in = obj.U(2, i);
                lf_in = car.lf;
                lr_in = car.lr;
                m_in = car.mass;
                mu_in = car.mu_s;
                omegaB_in = obj.States(6, i);
                varphi_in = obj.States(3, i);
                vx_in = obj.States(4, i);
                vy_in = obj.States(5, i);
                
                inputs = cell(15,1);
                inputs{i} = 
               
                if (abs(alphaf) < alphaFMax  && abs(alphar) < alphaRMax) % front grip, rear grip
                    obj.fn(:,i) = statesdot_fgrg(C_in, ...
                                                Flong_in, ...
                                                Fz_f_in, ...
                                                Fz_r_in, ...
                                                Iz_in, ...
                                                R_in, ...
                                                delta_in, ...
                                                lf_in, ...
                                                lr_in, ...
                                                m_in, ...
                                                mu_in, ...
                                                omegaB_in, ...
                                                varphi_in, ...
                                                vx_in, ...
                                                vy_in);

                    obj.A(:,i) = A_fgrg(C_in, ...
                                        Flong_in, ...
                                        Fz_f_in, ...
                                        Fz_r_in, ...
                                        Iz_in, ...
                                        R_in, ...
                                        delta_in, ...
                                        lf_in, ...
                                        lr_in, ...
                                        m_in, ...
                                        mu_in, ...
                                        omegaB_in, ...
                                        varphi_in, ...
                                        vx_in, ...
                                        vy_in);

                    obj.B(:,:,i) = B_fgrg(C_in, ...
                                        Flong_in, ...
                                        Fz_f_in, ...
                                        Fz_r_in, ...
                                        Iz_in, ...
                                        R_in, ...
                                        delta_in, ...
                                        lf_in, ...
                                        lr_in, ...
                                        m_in, ...
                                        mu_in, ...
                                        omegaB_in, ...
                                        varphi_in, ...
                                        vx_in, ...
                                        vy_in);
               
               
                elseif (abs(alphaf) >= alphaFMax && abs(alphar) < alphaRMax) % front slip, rear grip
                     obj.fn(:,i) = statesdot_fsrg(C_in, ...
                                                Flong_in, ...
                                                Fz_f_in, ...
                                                Fz_r_in, ...
                                                Iz_in, ...
                                                R_in, ...
                                                delta_in, ...
                                                lf_in, ...
                                                lr_in, ...
                                                m_in, ...
                                                mu_in, ...
                                                omegaB_in, ...
                                                varphi_in, ...
                                                vx_in, ...
                                                vy_in);

                    obj.A(:,i) = A_fsrg(C_in, ...
                                        Flong_in, ...
                                        Fz_f_in, ...
                                        Fz_r_in, ...
                                        Iz_in, ...
                                        R_in, ...
                                        delta_in, ...
                                        lf_in, ...
                                        lr_in, ...
                                        m_in, ...
                                        mu_in, ...
                                        omegaB_in, ...
                                        varphi_in, ...
                                        vx_in, ...
                                        vy_in);

                    obj.B(:,:,i) = B_fsrg(C_in, ...
                                        Flong_in, ...
                                        Fz_f_in, ...
                                        Fz_r_in, ...
                                        Iz_in, ...
                                        R_in, ...
                                        delta_in, ...
                                        lf_in, ...
                                        lr_in, ...
                                        m_in, ...
                                        mu_in, ...
                                        omegaB_in, ...
                                        varphi_in, ...
                                        vx_in, ...
                                        vy_in);
                   
                elseif (abs(alphaf) < alphaFMax && abs(alphar) >= alphaRMax) % front grip, rear slip
                     obj.fn(:,i) = statesdot_fgrs(C_in, ...
                                                Flong_in, ...
                                                Fz_f_in, ...
                                                Fz_r_in, ...
                                                Iz_in, ...
                                                R_in, ...
                                                delta_in, ...
                                                lf_in, ...
                                                lr_in, ...
                                                m_in, ...
                                                mu_in, ...
                                                omegaB_in, ...
                                                varphi_in, ...
                                                vx_in, ...
                                                vy_in);

                    obj.A(:,i) = A_fgrs(C_in, ...
                                        Flong_in, ...
                                        Fz_f_in, ...
                                        Fz_r_in, ...
                                        Iz_in, ...
                                        R_in, ...
                                        delta_in, ...
                                        lf_in, ...
                                        lr_in, ...
                                        m_in, ...
                                        mu_in, ...
                                        omegaB_in, ...
                                        varphi_in, ...
                                        vx_in, ...
                                        vy_in);

                    obj.B(:,:,i) = B_fgrs(C_in, ...
                                        Flong_in, ...
                                        Fz_f_in, ...
                                        Fz_r_in, ...
                                        Iz_in, ...
                                        R_in, ...
                                        delta_in, ...
                                        lf_in, ...
                                        lr_in, ...
                                        m_in, ...
                                        mu_in, ...
                                        omegaB_in, ...
                                        varphi_in, ...
                                        vx_in, ...
                                        vy_in);

                elseif (abs(alphaf) >= alphaFMax && abs(alphar) >= alphaRMax) % front slip, rear slip
                     obj.fn(:,i) = statesdot_fsrs(C_in, ...
                                                Flong_in, ...
                                                Fz_f_in, ...
                                                Fz_r_in, ...
                                                Iz_in, ...
                                                R_in, ...
                                                delta_in, ...
                                                lf_in, ...
                                                lr_in, ...
                                                m_in, ...
                                                mu_in, ...
                                                omegaB_in, ...
                                                varphi_in, ...
                                                vx_in, ...
                                                vy_in);

                    obj.A(:,i) = A_fsrs(C_in, ...
                                        Flong_in, ...
                                        Fz_f_in, ...
                                        Fz_r_in, ...
                                        Iz_in, ...
                                        R_in, ...
                                        delta_in, ...
                                        lf_in, ...
                                        lr_in, ...
                                        m_in, ...
                                        mu_in, ...
                                        omegaB_in, ...
                                        varphi_in, ...
                                        vx_in, ...
                                        vy_in);

                    obj.B(:,:,i) = B_fsrs(C_in, ...
                                        Flong_in, ...
                                        Fz_f_in, ...
                                        Fz_r_in, ...
                                        Iz_in, ...
                                        R_in, ...
                                        delta_in, ...
                                        lf_in, ...
                                        lr_in, ...
                                        m_in, ...
                                        mu_in, ...
                                        omegaB_in, ...
                                        varphi_in, ...
                                        vx_in, ...
                                        vy_in);
                end
           end
       end
       
       function cost = cost(obj,inVec)
           
       end
       
       
       
       
   end
       
    
end
