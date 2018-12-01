classdef MPC < handle
    properties
        
        A         % (6x6xN)
        B         % (6x2xN)
        
        % Optimization variables
        thetaA    % (1xN) 
        States    % (6xN) [X, Y, varphi, vx, vy, omegaB]'
        delta     % (1xN) 
        F_long    % (1xN) 
        v         % (1xN) 
        ddelta    % (1xN)
        dF_long   % (1xN)
        dv        % (1xN)
        
        nVarsPerIter % number of optimzation variables per time step     
        
        Ts        % sampling period
        N         % Horizon length
        
        slipflag
        
    end
    
    methods
        function obj = MPC(Ts,N, car)
            addpath ./TaylorEqs/;
            obj.Ts = Ts;
            obj.N = N;
            
            
            % set initial optimization variables
            obj.thetaA      = linspace(0,5,N);
            obj.States      = zeros(6,N); %[car.X, car.Y, car.phi, car.vx, car.vy, car.omegaB]';
            obj.States(1,:) = car.X;
            obj.States(2,:) = car.Y;
            obj.States(3,:) = car.phi;
            obj.States(4,:) = car.vx;
            obj.States(5,:) = car.vy;
            obj.States(6,:) = car.omegaB;
            obj.delta       = zeros(1,N); %[car.F_long, car.delta]';
            obj.F_long       = 400*ones(1,N); %[car.F_long, car.delta]';
            obj.v           = 10*ones(1,N);
            obj.ddelta      = zeros(1,N);
            obj.dF_long      = zeros(1,N);
            obj.dv          = zeros(1,N);
            
            obj.nVarsPerIter = 13;
            
            
        end
        
        function [] = setLinPoints(obj,States)
            N1 = obj.N;
            nVarsPerIter = obj.nVarsPerIter;
            
            U(:,1) = [obj.delta(N1); obj.F_long(N1)]; %[delta, F_long]
            
            % set initial states to current sim states
            obj.States(:,1) = States;
            
            for i = 2:N1-1
                obj.States(:,i) = obj.States(:,i+1);
            end

            obj.States(:,N1) = obj.A(:,:,N1)*obj.States(:,N1) + obj.B(:,:,N1)*U;


            % shift optimization vars by one sampling period to use as
            % linearization points and inital guess in fmincon
            for i = 1:N1-1
                obj.F_long(i) = obj.F_long(i+1);
                obj.delta(i) = obj.delta(i+1);
                obj.v(i) = obj.v(i+1);
                obj.thetaA(i) = obj.thetaA(i+1);
                obj.ddelta(i) = obj.ddelta(i+1);
                obj.dF_long(i) = obj.dF_long(i+1);
                obj.dv(i) = obj.dv(i+1);
            end
            
        end
        
        function [] = linearizeModel(obj, car)
            
            Fz_f = car.mass*9.81*car.dm;
            Fz_r = car.mass*9.81*(1-car.dm);
            
            N1 = obj.N;
            for i = 1:N1
                lf = car.dm*car.length;  %lengths from cg to front and rear end
                lr = (1-car.dm)*car.length;
                R = car.mu_k/car.mu_s;
                alphaFMax = atan2(3*car.mu_s*Fz_f, car.C);
                alphaRMax = atan2(3*car.mu_s*Fz_r, car.C);
                
                alphaf = atan2(obj.States(5, i) + lr*obj.States(6, i), obj.States(4, i)) - obj.delta(i);
                alphar = atan2(obj.States(5, i) - lf*obj.States(6, i), obj.States(4, i));
                
                %%%   INPUT VECTOR for f, A, B %%%%%
                %    C,F_long,Fz_f,Fz_r,Iz,R,delta,lf,lr,m,mu,omegaB,varphi,vx,vy
                C = car.C;
                F_long = obj.F_long(i);
                Fz_f = Fz_f;
                Fz_r = Fz_r;
                Iz = car.Iz;
                R = R;
                delta_in = obj.delta(i);
                lf = lf;
                lr = lr;
                m = car.mass;
                mu = car.mu_s;
                omegaB = obj.States(6, i);
                varphi = obj.States(3, i);
                vx = obj.States(4, i);
                vy = obj.States(5, i);
                
                if abs(vx) < 1 %10^-8
                    vx = 1; % 1e-8;
                end
                if abs(vy) < 10^-8
                    vy = 10^-8;
                end
                
                if (abs(alphaf) < alphaFMax  && abs(alphar) < alphaRMax) % front grip, rear grip
                    Ac = A_fgrg(C,Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                    Bc = B_fgrg(C,F_long,Fz_f,Iz,R,delta_in,lf,m,mu,omegaB,vx,vy);
                    obj.slipflag = 1;
                    
                elseif (abs(alphaf) >= alphaFMax && abs(alphar) < alphaRMax) % front slip, rear grip
                    Ac = A_fsrg(C,Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                    Bc = B_fsrg(F_long,Fz_f,Iz,R,delta_in,lf,m,mu,omegaB,vx,vy);
                    obj.slipflag = 2;
                    
                elseif (abs(alphaf) < alphaFMax && abs(alphar) >= alphaRMax) % front grip, rear slip
                    Ac = A_fgrs(C,Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                    Bc = B_fgrs(C,F_long,Fz_f,Iz,R,delta_in,lf,m,mu,omegaB,vx,vy);
                    obj.slipflag = 3;
                    
                elseif (abs(alphaf) >= alphaFMax && abs(alphar) >= alphaRMax) % front slip, rear slip
                    Ac = A_fsrs(Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                    Bc = B_fsrs(F_long,Fz_f,Iz,R,delta_in,lf,m,mu,omegaB,vx,vy);
                    obj.slipflag = 4;
                end
                
                % Convert To Discrete
                obj.A(:,:,i) = expm(Ac*obj.Ts);
                phitemp = zeros(size(Ac));
                for k = 0:1:10
                    phitemp = phitemp + Ac^k*obj.Ts^(k+1)/(factorial(k+1));
                end
                obj.B(:,:,i) = phitemp*Bc;

              
            end
        end
        
        function [G,H] = getCostMatrices(obj, ax, ay, bx, by, cx, cy, dx, dy)
            H = zeros(obj.nVarsPerIter*obj.N);
            G = zeros(obj.nVarsPerIter*obj.N,1);
            
            % set weights cost functions terms
            qc = 10;
            ql = 10;
            gamma = 1;
            rdelta = 100;
            rF_long = 1;
            rv = 1;
            
            for i = 1:obj.N
                X0 = obj.States(1,i);
                Y0 = obj.States(2,i);
                theta0 = obj.thetaA(i);
                ec0 = getEc(X0,Y0,ax,ay,bx,by,cx,cy,dx,dy,theta0);
                el0 = getEl(X0,Y0,ax,ay,bx,by,cx,cy,dx,dy,theta0);
                
                Gtemp = getDf(obj.Ts,X0,Y0,ax,ay,bx,by,cx,cy,dx,dy,ec0,el0,gamma,qc,ql,theta0);
                Htemp = getDDf(X0,Y0,ax,ay,bx,by,cx,cy,dx,dy,qc,ql,rF_long,rdelta,rv,theta0);
                
                G(1+(i-1)*obj.nVarsPerIter:i*obj.nVarsPerIter) = Gtemp;
                H(1+(i-1)*obj.nVarsPerIter:i*obj.nVarsPerIter,1+(i-1)*obj.nVarsPerIter:i*obj.nVarsPerIter) = Htemp;
                
            end
            
        end
        
        function cost = cost(obj,opt)
            % weights
            qc = 100;
            ql = 100;
            gamma = .001;
            ru1 = 1;
            ru2 = .1;
            rv = 1;
            R = [ru1, 0 , 0;
                0, ru2, 0;
                0, 0, rv];
            varsPerIter = 3;
            
            % inVec = [theta, F_long_1, delta_1, v_1, F_long_2, delta_2, v_2, ...
            lenInVec = length(opt);
            theta = opt(1);
            F_long_opt = opt(2:varsPerIter:lenInVec-2);
            delta_opt = opt(3:varsPerIter:lenInVec-1);
            v_opt = opt(4:varsPerIter:lenInVec);
            
            
            % these dont need to be outputed during optimization, they can
            % be computed and stored from the other opt variables after its all done
            thetaTemp = zeros(1,obj.N);
            StatesTemp = zeros(6,obj.N);
            thetaTemp(1) = theta;
            StatesTemp(:,1) = obj.States(:, 1);
            
            
            cost = -gamma*obj.Ts*v_opt(1);
            for k = 2:obj.N
                % step theta forward in time
                thetaTemp(k) = thetaTemp(k-1) + v_opt(k-1)*obj.Ts;
                
                % step the states forward in time
                StatesTemp(:,k) = StatesTemp(:,k-1)  ...
                    + obj.Ts*(obj.A(:,:,k-1)*(StatesTemp(:,k-1)-obj.States(:,k-1)) ...
                    + obj.B(:,:,k-1)*([F_long_opt(k-1),delta_opt(k-1)]' - [obj.F_long(k-1), obj.delta(k-1)]' ) + obj.f(:,k-1));
                
                % add the next term in the cost function
                cost = cost ...
                    + ql*[StatesTemp(1,k); StatesTemp(2,k); thetaTemp(k)]'*obj.Gamma_el(:,:,k)*[StatesTemp(1,k); StatesTemp(2,k); thetaTemp(k)] ...
                    + ql*obj.C_el(:,k)'*[StatesTemp(1,k); StatesTemp(2,k); thetaTemp(k)] ...
                    + qc*[StatesTemp(1,k); StatesTemp(2,k); thetaTemp(k)]'*obj.Gamma_ec(:,:,k)*[StatesTemp(1,k); StatesTemp(2,k); thetaTemp(k)] ...
                    + qc*obj.C_ec(:,k)'*[StatesTemp(1,k); StatesTemp(2,k); thetaTemp(k)] ...
                    - gamma*v_opt(k)*obj.Ts ...
                    + [F_long_opt(k)-F_long_opt(k-1); delta_opt(k)-delta_opt(k-1); v_opt(k)-v_opt(k-1)]'*R*[F_long_opt(k)-F_long_opt(k-1); delta_opt(k)-delta_opt(k-1); v_opt(k)-v_opt(k-1)];
                
            end
        end
        

        
        function [Aeq, Beq] = getEqualityCons(obj)
            nVarsPerIter = 13;
            nVars = nVarsPerIter*obj.N;
            
            Aeq = zeros(nVars);
            Beq = zeros(nVars, 1);
            % initialize theta and states
            Aeq(1:7, 1:7) = eye(7);
            Beq(1,1) = obj.thetaA(1);
            Beq(2:7) = obj.States(:,1);
            
            N1 = obj.N;
            
            r_start = 1 + nVarsPerIter;
            r_stop = 2*nVarsPerIter;
            c_start = 1;
            c_stop = 2*nVarsPerIter;
            
            
            for i = 2:N1
                temp = zeros(nVarsPerIter, 2*nVarsPerIter);
                % theta_k + vk*Ts - theta_k+1 = 0
                temp(1, 1) = 1;
                temp(1, 10) = obj.Ts;
                temp(1, nVarsPerIter+1) = -1;
                
                % A_k*X_k + B_k*U_k - X_k+1 = 0
                temp(2:7, 2:7) = obj.A(:,:,i-1);
                temp(2:7, 8:9) = obj.B(:,:,i-1);
                temp(2:7, 2+nVarsPerIter:7+nVarsPerIter) = -eye(size(obj.A(:,:,i)));
                
                % delta_k+1 - delta_k - ddelta_k+1 = 0
                temp(11, 8) = -1;
                temp(11, 8+nVarsPerIter) = 1;
                temp(11, 11+nVarsPerIter) = -1;
                
                % F_long_k+1 - F_long_k - dF_long_k+1 = 0
                temp(12, 9) = -1;
                temp(12, 9+nVarsPerIter) = 1;
                temp(12, 12+nVarsPerIter) = -1;
                
                % v_k+1 - v_k - dv_k+1 = 0
                temp(13, 10) = -1;
                temp(13, 10+nVarsPerIter) = 1;
                temp(13, 13+nVarsPerIter) = -1;
                

                % Put temp into big Aeq matrix
                Aeq(r_start:r_stop, c_start:c_stop) = temp;
                
                % re-calc indices to insert into for next iteration
                r_start = r_start + nVarsPerIter;
                r_stop = r_stop + nVarsPerIter;
                c_start = c_start + nVarsPerIter;
                c_stop = c_stop + nVarsPerIter;

            end
            
        end
        
        
        
        function [lb, ub] = getBounds(obj, car, rt)
            nn = obj.nVarsPerIter;
            nOptVars = obj.nVarsPerIter*obj.N;
            
            lb = zeros(nOptVars,1);
            ub = zeros(nOptVars,1);
            
            % theta constraint
            lb(1:nn:end) = 0;
            ub(1:nn:end) =  max(rt.theta);
            
            % X constraint
            lb(2:nn:end) = -100000;
            ub(2:nn:end) = 100000;
            
            % Y constraint
            lb(3:nn:end) = -100000;
            ub(3:nn:end) = 100000;
            
            % phi constraint
            lb(4:nn:end) = -100000;
            ub(4:nn:end) = 100000;
            % vx constraint
            lb(5:nn:end) = -100;
            ub(5:nn:end) = 100;
            
            % vy constraint
            lb(6:nn:end) = -100;
            ub(6:nn:end) = 100;
            
            % omegaB constraint
            lb(7:nn:end) = -100;
            ub(7:nn:end) = 100;
            
            % delta bounds
            lb(8:nn:end) = car.delta_min;
            ub(8:nn:end) = car.delta_max;
            
            % F_long bounds
            lb(9:nn:end) = car.F_long_brake;
            ub(9:nn:end) = car.F_long_accel;
            
            % v bounds
            lb(10:nn:end) = 0;
            ub(10:nn:end) = 100;
            
            % ddelta bounds
            lb(11:nn:end) = -100;
            ub(11:nn:end) = 100;
            
            % dF_long bounds
            lb(12:nn:end) = -1000;
            ub(12:nn:end) = 1000;
            
            % dv bounds
            lb(13:nn:end) = -100;
            ub(13:nn:end) = 100;
            
        end
        
    end
    
end
