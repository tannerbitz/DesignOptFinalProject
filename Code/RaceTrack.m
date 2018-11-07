classdef RaceTrack < handle
    properties
        controlPoints  % racetrack control points loaded from file
        progress       % race track progress varibale
        arcLen         % race track arc length corresponding to progress variable
        X              % x coordinate corresponding to progress variable
        Y              % y coordinate correspinding to progress variable
        width          % track width
        Xout           % outer edge x coordinate
        Yout           % outer edge y coordinate
        Xin            % inner edge x coordinate
        Yin            % inner edge y coordinate
        
    end
    
    methods
        
        function obj = RaceTrack(inputFile,width_in)
            obj.width = width_in;
            obj.controlPoints = importdata(inputFile);
            
        end
        
        function obj = computeRaceTrack(obj)
            order = 3; % cubic polynomial
            n =  size(obj.controlPoints,1); % number of control points
            
            % append on 4 points first 4 points, this make the track c0 and
            % c1 continuous at 1 point
            p = obj.controlPoints;
            p(n+1,:) = p(1,:);
            p(n+2,:) = p(2,:);
            p(n+3,:) = p(3,:);
            p(n+4,:) = p(4,:);
            
            n = size(p,1); % fix number of control points
            
            N = 10000;
            T = linspace(0,1,n-order+2);
            yt = linspace(1/(n-2),1-1/(n-2),N); % ignore first and last (1/n-2) values of yt, this starts and ends the track where it smoothly meets
            p_spl = obj.DEBOOR(T,p,yt,order);  % compute coordinates of b-spline track
            
            
            % compute length of each linear segment of track
            arcLen1 = zeros(size(p_spl,1),1);
            for i = 2:size(p_spl,1)
                arcLen1(i) = sqrt((p_spl(i-1,1)-p_spl(i,1))^2 + (p_spl(i-1,2)-p_spl(i,2))^2) + arcLen1(i-1);
            end
            
            
            % The sefments of the b-spline track are not equal arc length, so we
            % intperolate it onto a unifrom arc length vector
            N = 1000;
            arcLenMin = min(arcLen1);
            arcLenMax = max(arcLen1);
            obj.arcLen = linspace(arcLenMin,arcLenMax,N);
            obj.X = interp1(arcLen1,p_spl(:,1),obj.arcLen);
            obj.Y = interp1(arcLen1,p_spl(:,2),obj.arcLen);
            
            
            
            figure
            plot(obj.X,obj.Y, 'k--')
            hold on
            axis equal;
            
            % compute track edge locations
            for i = 1:N
                if i == 1
                    nVec(i,:) = [-(obj.Y(i+1)-obj.Y(N-2)), (obj.X(i+1)-obj.X(N-2))];
                    nVec(i,:) = nVec(i,:)/norm(nVec(i,:));
                elseif i == N
                    nVec(i,:) = [-(obj.Y(2)-obj.Y(i-1)), (obj.X(2)-obj.X(i-1))];
                    nVec(i,:) = nVec(i,:)/norm(nVec(i,:));
                else
                    nVec(i,:) = [-(obj.Y(i+1)-obj.Y(i-1)), (obj.X(i+1)-obj.X(i-1))];
                    nVec(i,:) = nVec(i,:)/norm(nVec(i,:));
                end
                
                obj.Xout(i) = obj.X(i) - nVec(i,1)*obj.width;
                obj.Yout(i) = obj.Y(i) - nVec(i,2)*obj.width;
                
                obj.Xin(i) = obj.X(i) + nVec(i,1)*obj.width;
                obj.Yin(i) = obj.Y(i) + nVec(i,2)*obj.width;
                
            end
            
            plot(obj.Xin,obj.Yin, 'r-', obj.Xout,obj.Yout ,'r-');
            
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