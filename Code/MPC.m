classdef MPC < handle
    properties
        f         % (6xN)
        A         % (6x6xN)
        B         % (6x2xN)
        Gamma_el  % (3x3xN)
        Gamma_ec  % (3x3xN)
        C_el      % (3xN)
        C_ec      % (3xN)
        
        States    % (6xN) [X, Y, varphi, vx, vy, omegaB]'
        thetaA    % (1xN) [thetaA]
        F_long    % (1xN) [Flong]
        delta     % (1xN) [delta]
        v         % (1xN) [v]
        
        rt
        car
        
        
        Ts        % sampling period
        N         % Horizon length
        
    end
    
    methods
        function obj = MPC(Ts,N, X0, Y0 ,rt,car)
            addpath ./TaylorEqs/;
            obj.Ts = Ts;
            obj.N = N;
            
            
            % set initial predicition horizon
            obj.States =  zeros(6,N); %[car.X, car.Y, car.phi, car.vx, car.vy, car.omegaB]';
            obj.States(1,:) = X0;
            obj.States(2,:) = Y0;
            obj.States(4,:) = 10;
            obj.F_long      = 400*ones(1,N); %[car.F_long, car.delta]';
            obj.delta       = zeros(1,N); %[car.F_long, car.delta]';
            obj.thetaA      = linspace(0,5,N);
            obj.v           = 10*ones(1,N);
            
            obj.rt = rt;
            obj.car = car;
            
            
        end
        
        function [] = setLinPoints(obj,States)
            N1 = obj.N;
            
            % set initial states to current sim states
            obj.States(:,1) = States;
            
            % shift states by one sampling period
            for i = 2:N1-1
                obj.States(:,i) = obj.States(:,i+1);
            end
            
            % compute last state
            obj.States(:,N1) = obj.States(:,N1) + obj.Ts*obj.f(:,N1);
            
            % shift optimization vars by one sampling period to use as
            % linearization points and inital guess in fmincon
            for i = 1:N1-1
                obj.F_long(i) = obj.F_long(i+1);
                obj.delta(i) = obj.delta(i+1);
                obj.v(i) = obj.v(i+1);
                obj.thetaA(i) = obj.thetaA(i+1);
            end
            
        end
        
        function [] = linearize(obj, car, ax, ay, bx, by, cx, cy, dx, dy)
            
            Fz_f = car.mass*9.81*car.dm;
            Fz_r = car.mass*9.81*(1-car.dm);
            
            N1 = obj.N;
            for i = 1:N1
                %%%   INPUT VECTOR for Gamma_ec, Gamma_el, Cec, Cel %%%%%
                %    (X,Y,ax,ay,bx,by,cx,cy,dx,dy,theta)
                X_in = obj.States(1, i);
                Y_in = obj.States(2, i);
                thetaA_in = obj.thetaA(i);
                
                inputs = {X_in, Y_in, ax, ay, bx, by, cx, cy, dx, dy, thetaA_in};
                obj.Gamma_ec(:, :, i) = GammaEc(inputs{:});
                obj.Gamma_el(:, :, i) = GammaEl(inputs{:});
                obj.C_ec(:, i) = Cec(inputs{:});
                obj.C_el(:, i) = Cel(inputs{:});
                
                lf = car.dm*car.length;  %lengths from cg to front and rear end
                lr = (1-car.dm)*car.length;
                R = car.mu_k/car.mu_s;
                alphaFMax = atan2(3*car.mu_s*Fz_f, car.C);
                alphaRMax = atan2(3*car.mu_s*Fz_r, car.C);
                
                alphaf = atan2(obj.States(5, i) + lr*obj.States(6, i), obj.States(4, i)) - obj.delta(i);
                alphar = atan2(obj.States(5, i) - lf*obj.States(6, i), obj.States(4, i));
                
                %%%   INPUT VECTOR for f, A, B %%%%%
                %    C,Flong,Fz_f,Fz_r,Iz,R,delta,lf,lr,m,mu,omegaB,varphi,vx,vy
                C = car.C;
                Flong = obj.F_long(i);
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
                    obj.f(:,i) = statesdot_fgrg(C,Flong,Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                    obj.A(:,:,i) = A_fgrg(C,Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                    obj.B(:,:,i) = B_fgrg(C,Flong,Fz_f,Iz,R,delta_in,lf,lr,m,mu,omegaB,vx,vy);
                    
                elseif (abs(alphaf) >= alphaFMax && abs(alphar) < alphaRMax) % front slip, rear grip
                    obj.f(:,i) = statesdot_fsrg(C,Flong,Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                    obj.A(:,:,i) = A_fsrg(C,Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                    obj.B(:,:,i) = B_fsrg(Flong,Fz_f,Iz,delta_in,lf,lr,m,mu,omegaB,vx,vy);
                    
                elseif (abs(alphaf) < alphaFMax && abs(alphar) >= alphaRMax) % front grip, rear slip
                    obj.f(:,i) = statesdot_fgrs(C,Flong,Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                    obj.A(:,:,i) = A_fgrs(C,Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                    obj.B(:,:,i) = B_fgrs(C,Flong,Fz_f,Iz,R,delta_in,lf,lr,m,mu,omegaB,vx,vy);
                    
                elseif (abs(alphaf) >= alphaFMax && abs(alphar) >= alphaRMax) % front slip, rear slip
                    obj.f(:,i) = statesdot_fsrs(Flong,Fz_f,Fz_r,Iz,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                    obj.A(:,:,i) = A_fsrs(Fz_f,Fz_r,Iz,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                    obj.B(:,:,i) = B_fsrs(Flong,Fz_f,Iz,delta_in,lf,lr,m,mu,omegaB,vx,vy);
                end
            end
        end
        
        function cost = cost(obj,opt)
            % weights
            qc = 1;
            ql = 1;
            gamma = 1;
            ru1 = 1;
            ru2 = 1;
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
        
        function cost = cost2(obj,opt)
            % weights
            qc = 100;
            ql = 10;
            gamma = 1;
            ru1 = 1;
            ru2 = 1;
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
                    + obj.Ts*getRHS(obj,StatesTemp(:,k-1), F_long_opt(k-1), delta_opt(k-1));
                
                % compute errors
                [ec,el] = obj.rt.getErrors(StatesTemp(1,k),StatesTemp(2,k),thetaTemp(k));
                
                
                % add the next term in the cost function
                cost = cost ...
                    + ql*el^2 + qc*ec^2  ...
                    - gamma*v_opt(k)*obj.Ts ...
                    + [F_long_opt(k)-F_long_opt(k-1); delta_opt(k)-delta_opt(k-1); v_opt(k)-v_opt(k-1)]'*R*[F_long_opt(k)-F_long_opt(k-1); delta_opt(k)-delta_opt(k-1); v_opt(k)-v_opt(k-1)];
                
            end
        end
        
        function States_dot = getRHS(obj,States_in, F_long_in, delta_in)
            Fz_f = obj.car.mass*9.81*obj.car.dm;
            Fz_r = obj.car.mass*9.81*(1-obj.car.dm);
            
            
            lf = obj.car.dm*obj.car.length;  %lengths from cg to front and rear end
            lr = (1-obj.car.dm)*obj.car.length;
            R = obj.car.mu_k/obj.car.mu_s;
            alphaFMax = atan2(3*obj.car.mu_s*Fz_f, obj.car.C);
            alphaRMax = atan2(3*obj.car.mu_s*Fz_r, obj.car.C);
            
            alphaf = atan2(States_in(5) + lr*States_in(6), States_in(4)) - delta_in;
            alphar = atan2(States_in(5) - lf*States_in(6), States_in(4));
            
            
            %%%   INPUT VECTOR for f, A, B %%%%%
            %    C,Flong,Fz_f,Fz_r,Iz,R,delta,lf,lr,m,mu,omegaB,varphi,vx,vy
            C = obj.car.C;
            Flong = F_long_in;
            Fz_f = Fz_f;
            Fz_r = Fz_r;
            Iz = obj.car.Iz;
            R = R;
            delta_in = delta_in;
            lf = lf;
            lr = lr;
            m = obj.car.mass;
            mu = obj.car.mu_s;
            omegaB = States_in(6);
            varphi = States_in(3);
            vx = States_in(4);
            vy = States_in(5);
            
            if abs(vx) < 1 %10^-8
                vx = 1; % 1e-8;
            end
            if abs(vy) < 10^-8
                vy = 10^-8;
            end
            
            if (abs(alphaf) < alphaFMax  && abs(alphar) < alphaRMax) % front grip, rear grip
                States_dot = statesdot_fgrg(C,Flong,Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                
            elseif (abs(alphaf) >= alphaFMax && abs(alphar) < alphaRMax) % front slip, rear grip
                States_dot = statesdot_fsrg(C,Flong,Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                
            elseif (abs(alphaf) < alphaFMax && abs(alphar) >= alphaRMax) % front grip, rear slip
                States_dot = statesdot_fgrs(C,Flong,Fz_f,Fz_r,Iz,R,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
                
            elseif (abs(alphaf) >= alphaFMax && abs(alphar) >= alphaRMax) % front slip, rear slip
                States_dot = statesdot_fsrs(Flong,Fz_f,Fz_r,Iz,delta_in,lf,lr,m,mu,omegaB,varphi,vx,vy);
            end
        end
       
    
        function [] = setStates(obj,opt)
            varsPerIter = 3;
            % inVec = [theta, F_long_1, delta_1, v_1, F_long_2, delta_2, v_2, ...
            lenInVec = length(opt);
            theta = opt(1);
            F_long_opt = opt(2:varsPerIter:lenInVec-2);
            delta_opt = opt(3:varsPerIter:lenInVec-1);
            v_opt = opt(4:varsPerIter:lenInVec);
            
            thetaTemp = zeros(1,obj.N);
            StatesTemp = zeros(6,obj.N);
            thetaTemp(1) = theta;
            StatesTemp(:,1) = obj.States(1);
            for k = 2:obj.N
                % step theta forward in time
                thetaTemp(k) = thetaTemp(k-1) + v_opt(k)*obj.Ts;
                
                % step the states forward in time
                StatesTemp(:,k) = StatesTemp(:,k-1)  ...
                    + obj.Ts*(obj.A(:,:,k-1)*(StatesTemp(:,k-1)-obj.States(:,k-1)) ...
                    + obj.B(:,:,k-1)*([F_long_opt(k-1),delta_opt(k-1)]' - [obj.F_long(k-1)', obj.delta(k-1)]' ) + obj.f(:,k-1));
            end
            
            obj.States = StatesTemp;
            obj.thetaA = thetaTemp;
            obj.F_long = F_long_opt;
            obj.delta  = delta_opt;
            obj.v      = v_opt;
            
        end
        
        
        
        function [Aineq, Bineq] = getIneqCons(obj, car)
            N1 = obj.N;
            nVarsPerIter = 3;
            nOptVars = 1+nVarsPerIter*N1;
            Aineq = [];
            Bineq = [];
            
            % theta constraint
            temp1 = zeros(1, nOptVars);
            temp1(1) = -1;
            Aineq = [Aineq; temp1];
            Bineq = [Bineq; 0];
            
            % Flong constraints
            for i = 2:nVarsPerIter:nOptVars
                temp1 = zeros(1, nOptVars);
                temp2 = zeros(1, nOptVars);
                temp1(i) = 1;
                temp2(i) = -1;
                Aineq = [Aineq; temp1; temp2];
                Bineq = [Bineq; car.F_long_accel; -car.F_long_brake];
            end
            
            % delta constraints
            for i = 3:nVarsPerIter:nOptVars
                temp1 = zeros(1, nOptVars);
                temp2 = zeros(1, nOptVars);
                temp1(i) = 1;
                temp2(i) = -1;
                Aineq = [Aineq; temp1; temp2];
                Bineq = [Bineq; car.delta_max; -car.delta_min];
            end
            
            % v constraints
            for i = 4:nVarsPerIter:nOptVars
                temp1 = zeros(1, nOptVars);
                temp2 = zeros(1, nOptVars);
                temp1(i) = 1;
                temp2(i) = -1;
                Aineq = [Aineq; temp1; temp2];
                Bineq = [Bineq; car.v_max; -car.v_min];
            end
            
        end
        
        
        function [lb, ub] = getBounds(obj, car, rt)
            N1 = obj.N;
            nn = 3;
            nOptVars = 1+N1*nn;
            
            lb = zeros(nOptVars,1);
            ub = zeros(nOptVars,1);
            
            % theta constraint
            lb(1) = 0;
            ub(1) = max(rt.theta);
            
            % F_long bounds
            lb(2:nn:end) = car.F_long_brake;
            ub(2:nn:end) = car.F_long_accel;
            
            % delta bounds
            lb(3:nn:end) = car.delta_min;
            ub(3:nn:end) = car.delta_max;
            
            % v bounds
            lb(4:nn:end) = 0;
            ub(4:nn:end) = 100;
            
        end
        
    end
    
end
