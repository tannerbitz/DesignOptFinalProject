classdef Fiala
    %FIALA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radiusWheel    % wheel radius
        C           % coefficient of friction (Cr = Cf = C)
    end
    
    methods
        function obj = Fiala(radiusWheel_, C_)
            %FIALA Construct an instance of this class
            %   Detailed explanation goes here
            obj.radiusWheel = radiusWheel_;
            obj.C = C_;
        end
        
        function  fLong = getFLong(obj)
            %getFLong calculates the longitudinal foce on the tire
            fLong = something;
        end
        
        function fLat = getFLat(obj)
            %getFLat calculates the lateral force on the a tire
            fLat = something;
        end
    end
end

