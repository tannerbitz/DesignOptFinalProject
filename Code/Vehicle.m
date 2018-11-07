classdef Vehicle < handle
    %This provides a base class for all vehicles
    
    properties
        massVehicle     % vehicle mass
        vx              % velocity of vehicle in the inertial x direction
                        % (which is the body's longitudinal direction)
        vy              % velocity of vehicle in the inertial y direction 
                        % (which is the body's lateral direction)
        phi             % the angle between the global x direction and the
                        % inertial x direction (measured counter-clockwise)
        x               % global x coordinate
        y               % global y coordinate
        omegaB          % angular velocity of the vehicle body
    end
    
    methods
        function obj = Vehicle()
            %Constructor for Vehicle class
        end
        
        function obj = setupvehicle(obj, varargin)
            % Vehicle property defaults
            defaultMassVehicle = 0; %kg
            defaultVx = 0;          %m/s
            defaultVy = 0;          %m/s
            defaultPhi = 0;         %rads
            defaultX = 0;           %m
            defaultY = 0;           %m
            defaultOmegaB = 0;      %rads/s
            
            % Set Parser/Parse Inputs
            p = inputParser;
            addParameter(p, 'massVehicle', defaultMassVehicle);
            addParameter(p, 'vx', defaultVx);
            addParameter(p, 'vy', defaultVy);
            addParameter(p, 'phi', defaultPhi);
            addParameter(p, 'x', defaultX);
            addParameter(p, 'y', defaultY);
            addParameter(p, 'omegaB', defaultOmegaB);
            parse(p, varargin{:});
            
            % Set object properties
            obj.massVehicle = p.Results.massVehicle;
            obj.vx = p.Results.vx;
            obj.vy = p.Results.vy;
            obj.phi = p.Results.phi;
            obj.x = p.Results.x;
            obj.y = p.Results.y;
            obj.omegaB = p.Results.omegaB;
        end
        
               
    end
end

