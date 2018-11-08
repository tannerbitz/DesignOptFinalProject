classdef Fiala < handle
    %FIALA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radiusWheel    % wheel radius
        C              % tire stiffness (Cr = Cf = C)
        mu             % coefficient of friction
    end
    
    methods
        function obj = Fiala()
            %FIALA Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function obj = setupfiala(obj, varargin)
            % Fiala property defaults
            defaultRadiusWheel = 1; %m
            defaultC = 1;           %unitless
            
            % Set Parser/Parse Inputs
            p = inputParser;
            addParameter(p, 'radiusWheel', defaultRadiusWheel);
            addParameter(p, 'C', defaultC);
            addParameter(p, 'mu', defaultMu);
            parse(p, varargin{:});
            
            % Set object properties
            obj.massVehicle = p.Results.radiusWheel;
            obj.vx = p.Results.C;
            obj.mu = p.Results.mu;
        end
        
        function  [fLong, fLat] = getTireForces(obj, vehicleObj, Fz, accelBrakeFlag)
            %getTireForces calculates and returns the latitudinal and
            %longitudinal forces on a wheel. The inputs are the normal
            %force on the tire, Fz, and a flag that specifies whether the
            %car is accelerating or braking. 
            
            % Parse Inputs
            p = inputParser;
            addRequired(p, 'vehicleObj');
            addRequired(p, 'Fz');
            addRequired(p, 'accelBrakeFlag');
            parse(p, vehicleObj, Fz, accelBrakeFlag);
            
            vObjTemp = p.Results.vehicleObj;
            fzTemp = p.Results.Fz;
            accelBrakeFlag = p.Results.accelBrakeFlag;
            
            % Calculate sigmaLat/ sigmaLong/ sigma
            if strcmp(accelBrakeFlag, 'accel')
                sigmaLong = (obj.radiusWheel*vObjTemp.omegaW - vObjTemp.vx)/(obj.radiusWheel*vObjTemp.omegaW);
            elseif strcmp(accelBrakeFlag, 'brake')
                sigmaLong = (obj.radiusWheel*vObjTemp.omegaW - vObjTemp.vx)/(vObjTemp.vx);
            end
            
            
            sigmaLat = vObjTemp.vx/(obj.radiusWheel*vObjTemp.omegaW)*tan(vObjTemp.alpha);
            sigma = sqrt(sigmaLat^2 + sigmaLong^2);
            
            % Calculate theta
            theta = obj.C/(3*obj.mu*fzTemp);
            
            % Determine if tire is sliding/ total force
            thetaPoly = 3*theta*sigma - 3*(theta*sigma)^2 + (theta*sigma)^3;
            if thetaPoly < 1 %not sliding
                Ft = obj.mu*fzTemp*thetaPoly;
            else
                Ft = obj.mu*fzTemp;
            end
            
            % Calculate lat and long tire forces
            fLong = sigmaLong/sigma*Ft;
            fLat = sigmaLat/sigma*Ft;
            
        end
    end
end

