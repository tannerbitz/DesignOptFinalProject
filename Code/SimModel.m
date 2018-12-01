classdef SimModel < handle
    properties
        % car
        vehicle

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


        function updateOneTrack(obj,deltaT, delta, F_long)
            % update vehicle states using ODE45

            States0 = [obj.X, obj.Y, obj.phi, obj.vx, obj.vy, obj.omegaB]';


            [~,newStates] = ode45(@(t,States) obj.calcRhsOneTrack(t,States, delta,F_long), [0 deltaT], States0);

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
            Ffl_z = .5*obj.mass*9.81*(1-obj.dm) - obj.Fx*obj.hcg/obj.length - obj.Fy*obj.hcg/obj.width;
            Ffr_z = .5*obj.mass*9.81*(1-obj.dm) - obj.Fx*obj.hcg/obj.length + obj.Fy*obj.hcg/obj.width;
            Frl_z = .5*obj.mass*9.81*(obj.dm) + obj.Fx*obj.hcg/obj.length - obj.Fy*obj.hcg/obj.width;
            Frr_z = .5*obj.mass*9.81*(obj.dm) + obj.Fx*obj.hcg/obj.length + obj.Fy*obj.hcg/obj.width;

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
            omegaB_dot = 1/obj.Iz*( (Ffl_long+Ffr_long)*sin(delta)*lf + (Ffr_long-Ffl_long)*cos(delta)*obj.width/2 ...
                                   +(Ffl_lat+Ffr_lat)*cos(delta)*lf + (Ffl_lat-Ffr_lat)*sin(delta)*obj.width/2 ...
                                   - (Frr_lat + Frl_lat)*lr + (Frr_long+Frl_long)*obj.width/2 );

            States_dot = [X_dot,Y_dot, phi_dot, vx_dot, vy_dot, omegaB_dot]';

        end


        function [States_dot] = calcRhsOneTrack(obj,~,States, delta, F_long)

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
            Ff_z = obj.mass*9.81*(1 - obj.dm);
            Fr_z = obj.mass*9.81*obj.dm;

            % compute tire longitudenal and lateral forces with Fiala tire
            % model
            Ff_lat = obj.getTireForcesMassera(omegaB_in, delta, vy_in, vx_in, Ff_z, lf);
            Fr_lat = obj.getTireForcesMassera(omegaB_in, 0, vy_in, vx_in, Fr_z, -lr);

            Ff_long = F_long;
            Fr_long = F_long;

            % update Fx and Fy so it can be used in the next timestep;
            obj.Fx = Fr_long + - Ff_lat*sin(delta) + Ff_long*cos(delta);
            obj.Fy = Fr_lat + Ff_lat*cos(delta) + Ff_long*sin(delta);

            % compute derivatives of vehicle states
            X_dot = vx_in*cos(phi_in) - vy_in*sin(phi_in);
            Y_dot = vx_in*sin(phi_in) + vy_in*cos(phi_in);
            phi_dot = omegaB_in;
            vx_dot  = 1/obj.mass*(obj.Fx + obj.mass*vy_in*obj.omegaB);
            vy_dot  = 1/obj.mass*(obj.Fy - obj.mass*vx_in*obj.omegaB);
            omegaB_dot = 1/obj.Iz*(Ff_lat*lf*cos(delta) + Ff_long*lf*sin(delta) - Fr_lat*lr);


            States_dot = [X_dot,Y_dot, phi_dot, vx_dot, vy_dot, omegaB_dot]';

        end


        function [] = updateOneTrackLinear(obj, deltaT, delta, F_long)

            addpath ./TaylorEqs/

            Fz_f = obj.mass*9.81*obj.dm;
            Fz_r = obj.mass*9.81*(1-obj.dm);

            R = obj.mu_k/obj.mu_s;
            alphaFMax = atan2(3*obj.mu_s*Fz_f, obj.C);
            alphaRMax = atan2(3*obj.mu_s*Fz_r, obj.C);

            % compute lengths from cg to front and rear end
            lf = obj.dm*obj.length;
            lr = (1-obj.dm)*obj.length;

            alphaf = atan2(obj.vy + lr*obj.omegaB, obj.vx) - delta;
            alphar = atan2(obj.vy - lf*obj.omegaB, obj.vx);

            vx = obj.vx;
            vy = obj.vy;
            if abs(obj.vx) < 1 %10^-8
                vx = 1; % 1e-8;
            end
            if abs(obj.vy) < 10^-8
                vy = 10^-8;
            end


            if (abs(alphaf) < alphaFMax  && abs(alphar) < alphaRMax) % front grip, rear grip
                f = statesdot_fgrg(C,Flong,Fz_f,Fz_r,Iz,R,delta,lf,lr,m,mu,omegaB,varphi,vx,vy);
                Ac = A_fgrg(obj.C,Fz_f,Fz_r,obj.Iz,R,delta,lf,lr,obj.mass,obj.mu_s,obj.omegaB,obj.phi,vx,vy);
                               %C,Fz_f,Fz_r,    Iz,R,delta,lf,lr,       m,      mu,    omegaB, varphi,vx,vy)
                Bc = B_fgrg(obj.C,F_long,Fz_f,obj.Iz,R,delta,lf,obj.mass,obj.mu_s,obj.omegaB,vx,vy);
                               %C, Flong,Fz_f,    Iz,R,delta,lf,       m,    mu,      omegaB,vx,vy

            elseif (abs(alphaf) >= alphaFMax && abs(alphar) < alphaRMax) % front slip, rear grip
                Ac = A_fsrg(obj.C,Fz_f,Fz_r,obj.Iz,R,delta,lf,lr,obj.mass,obj.mu_s,obj.omegaB,obj.phi,vx,vy);
                               %C,Fz_f,Fz_r,    Iz,R,delta,lf,lr,       m,    mu,      omegaB, varphi,vx,vy
                Bc = B_fsrg(F_long,Fz_f,obj.Iz,R,delta,lf,obj.mass,obj.mu_s,obj.omegaB,vx,vy);
                           % Flong,Fz_f,    Iz,R,delta,lf,       m,      mu,    omegaB,vx,vy

            elseif (abs(alphaf) < alphaFMax && abs(alphar) >= alphaRMax) % front grip, rear slip
                Ac = A_fgrs(obj.C,Fz_f,Fz_r,obj.Iz,R,delta,lf,lr,obj.mass,obj.mu_s,obj.omegaB,obj.phi,vx,vy);
                            %   C,Fz_f,Fz_r,    Iz,R,delta,lf,lr,       m,    mu,      omegaB, varphi,vx,vy
                Bc = B_fgrs(obj.C,F_long,Fz_f,obj.Iz,R,delta,lf,obj.mass,obj.mu_s,obj.omegaB,vx,vy);
                             %  C, Flong,Fz_f,    Iz,R,delta,lf,       m,    mu,      omegaB,vx,vy

            elseif (abs(alphaf) >= alphaFMax && abs(alphar) >= alphaRMax) % front slip, rear slip
                Ac = A_fsrs(Fz_f,Fz_r,obj.Iz,R,delta,lf,lr,obj.mass,obj.mu_s,obj.omegaB,obj.phi,vx,vy);
                           %Fz_f,Fz_r,    Iz,R,delta,lf,lr,       m,      mu,    omegaB, varphi,vx,vy
                Bc = B_fsrs(F_long,Fz_f,obj.Iz,R,delta,lf,obj.mass,obj.mu_s,obj.omegaB,vx,vy);
                            %Flong,Fz_f,    Iz,R,delta,lf,       m,    mu,      omegaB,vx,vy

            end

            % Convert To Discrete
            states = [obj.X; obj.Y; obj.phi; obj.vx; obj.vy; obj.omegaB];
            inputs = [delta; F_long];

            %             A = expm(Ac*deltaT);
            %             phitemp = zeros(size(Ac));
            %             for k = 0:1:10
            %                 phitemp = phitemp + Ac^k*deltaT^(k+1)/(factorial(k+1));
            %             end
            %             B = phitemp*Bc;
            %
            %             states = A*states + B*inputs;

            states = states + deltaT*(Ac*states + Bc*inputs);

            obj.X = states(1);
            obj.Y = states(2);
            obj.phi = states(3);
            obj.vx = states(4);
            obj.vy = states(5);
            obj.omegaB = states(6);
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
                F_lat = -f_alpha + (2-R_mu)/(3*obj.mu_s*Fz_in)*abs(f_alpha)*f_alpha ...
                    - (1-2/3*R_mu)/(3*obj.mu_s*Fz_in)^2*f_alpha^3;
            else
                F_lat = -sign(alpha)*obj.mu_s*R_mu*Fz_in;
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
