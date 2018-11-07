classdef RaceTrack
    properties
        controlPoints  % racetrack control points to load from file
        progress       % race track progress varibale
        coordinates    % coordiates [x,y] corresponding to progress variable
    end
    
    methods
        
        function obj = RaceTrack(inputFile)
            obj.controlPoints = importdata(inputFile);
            obj.computeRaceTrack();
        end
        
        
        function obj = computeRaceTrack(obj)
            order = 3;
            n =  size(obj.controlPoints,1);
            
            p = obj.controlPoints;
            p(n+1,:) = p(1,:);
            p(n+2,:) = p(2,:);
            p(n+3,:) = p(3,:);
            p(n+4,:) = p(4,:);
            
            n = size(p,1);
            
            T = linspace(0,1,n-order+2);
            y = linspace(1/(n-2),1-1/(n-2),100000);
            p_spl = obj.DEBOOR(T,p,y,order);
            plot(p_spl(:,1),p_spl(:,2),'b.-','LineWidth',1);
            
            
        end
        
        function val = DEBOOR(obj,T,p,y,order)
            
            % function val = DEBOOR(T,p,y,order)
            %
            % INPUT:  T     St�tzstellen
            %         p     Kontrollpunkte (nx2-Matrix)
            %         y     Auswertungspunkte (Spaltenvektor)
            %         order Spline-Ordnung
            %
            % OUTPUT: val   Werte des B-Splines an y (mx2-Matrix)
            %
            % Date:   2007-11-27
            % Author: Jonas Ballani
            p;
            m = size(p,1);
            n = length(y);
            X = zeros(order,order);
            Y = zeros(order,order);
            a = T(1);
            b = T(end);
            T = [ones(1,order-1)*a,T,ones(1,order-1)*b];
            
            
            for l = 1:n
                t0 = y(l);
                id = find(t0 >= T);
                k = id(end);
                if (k > m)
                    return;
                end
                X(:,1) = p(k-order+1:k,1);
                Y(:,1) = p(k-order+1:k,2);
                
                for i = 2:order
                    for j = i:order
                        num = t0-T(k-order+j);
                        if num == 0
                            weight = 0;
                        else
                            s = T(k+j-i+1)-T(k-order+j);
                            weight = num/s;
                        end
                        X(j,i) = (1-weight)*X(j-1,i-1) + weight*X(j,i-1);
                        Y(j,i) = (1-weight)*Y(j-1,i-1) + weight*Y(j,i-1);
                    end
                end
                val(l,1) = X(order,order);
                val(l,2) = Y(order,order);
            end
            
        end
        
    end
    
end