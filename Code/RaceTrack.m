classdef RaceTrack < handle
    properties
        controlPoints  % racetrack control points loaded from file
        N1
        N2
        theta          % arc lengths along track at corresponding XY
        X              % x coordinates
        Y              % y coordinates
        width          % track width
        Xout           % outer edge x coordinate
        Yout           % outer edge y coordinate
        Xin            % inner edge x coordinate
        Yin            % inner edge y coordinate
    end

    methods

        function obj = RaceTrack(inputFile,width_in)
            % RACETRACK is the racetrack constructor. It takes in the width
            % of the track and an input file with x, y coordinators of 
            % control points to create the B-splines that make up the track
            obj.width = width_in;
            obj.controlPoints = importdata(inputFile);

        end
        
        function [theta_out,X_out,Y_out,lapCount] = getXY(obj,theta_in, lapCount)
            
            % if we're passed the end of the track, start at the beginning
            % again
            
            if theta_in(1) >= obj.theta(obj.N1)
                theta_out = theta_in - obj.theta(end)/2;
                lapCount = lapCount + 1;
            else
                theta_out = theta_in;
            end
            
            iStart = find(obj.theta <= theta_out(1));
            iStart = iStart(end);
            iEnd = find(obj.theta >= theta_out(end));
            iEnd = iEnd(1);
            
            X_out = interp1(obj.theta(iStart:iEnd), obj.X(iStart:iEnd),theta_out);
            Y_out = interp1(obj.theta(iStart:iEnd), obj.Y(iStart:iEnd),theta_out);
        end
        
        function [ax,bx,cx,dx,ay,by,cy,dy] = getCubicPolynomial(obj,theta_in,X_in,Y_in)
            
            % compute vacndermond matrix
            for i = 1:length(theta_in)
                V(i,:) = [1 theta_in(i) theta_in(i)^2 theta_in(i)^3];
            end
            
            
            coef_x = inv(V'*V)*V'*X_in';
            ax = coef_x(1);
            bx = coef_x(2);
            cx = coef_x(3);
            dx = coef_x(4);
            
            coef_y = inv(V'*V)*V'*Y_in';
            ay = coef_y(1);
            by = coef_y(2);
            cy = coef_y(3);
            dy = coef_y(4);
            
        end
        
        function obj = computeRaceTrack(obj)
            order = 3; % cubic polynomial
            n =  size(obj.controlPoints,1); % number of control points

            % append last 4 points as first 4 points, this make the track c0 and
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
            theta1 = zeros(size(p_spl,1),1);
            for i = 2:size(p_spl,1)
                theta1(i) = sqrt((p_spl(i-1,1)-p_spl(i,1))^2 + (p_spl(i-1,2)-p_spl(i,2))^2) + theta1(i-1);
            end


            % The sefments of the b-spline track are not equal arc length, so we
            % intperolate it onto a unifrom arc length vector
            
            obj.N1 = 1000;
            thetaMin = min(theta1);
            thetaMax = max(theta1);
            obj.theta = linspace(thetaMin,thetaMax,obj.N1);
            obj.X = interp1(theta1,p_spl(:,1),obj.theta);
            obj.Y = interp1(theta1,p_spl(:,2),obj.theta);

            % compute track edge locations
            for i = 1:obj.N1
                if i == 1
                    nVec(i,:) = [-(obj.Y(i+1)-obj.Y(obj.N1-2)), (obj.X(i+1)-obj.X(obj.N1-2))];
                    nVec(i,:) = nVec(i,:)/norm(nVec(i,:));
                elseif i == obj.N1
                    nVec(i,:) = [-(obj.Y(2)-obj.Y(i-1)), (obj.X(2)-obj.X(i-1))];
                    nVec(i,:) = nVec(i,:)/norm(nVec(i,:));
                else
                    nVec(i,:) = [-(obj.Y(i+1)-obj.Y(i-1)), (obj.X(i+1)-obj.X(i-1))];
                    nVec(i,:) = nVec(i,:)/norm(nVec(i,:));
                end

                obj.Xout(i) = obj.X(i) - nVec(i,1)*obj.width/2;
                obj.Yout(i) = obj.Y(i) - nVec(i,2)*obj.width/2;

                obj.Xin(i) = obj.X(i) + nVec(i,1)*obj.width/2;
                obj.Yin(i) = obj.Y(i) + nVec(i,2)*obj.width/2;

            end
            
            
            % lets include a second loop around the track, this makes it
            % easier to to compute the cubic spline when the car gets to
            % the end of the track
            obj.theta = [obj.theta, obj.theta(2:end)+obj.theta(end)];
            obj.X = [obj.X, obj.X(2:end)];
            obj.Y = [obj.Y, obj.Y(2:end)];
            
            obj.N2 = length(obj.theta);
            
            plot(obj.X,obj.Y, 'k--')
            hold on
            plot(obj.Xin,obj.Yin, 'r-', obj.Xout,obj.Yout ,'r-');
%             plot(p(:,1),p(:,2),'bo-')
            
            axis equal;
            
        end
        
        function [ec, el, XA, YA] = getErrors(obj,X,Y,thetaA)
            if thetaA < 10^-8
                i1 = 1;
            else
                i1 = find(obj.theta <= thetaA);
                i1 = i1(end);
            end
            i2 = find(obj.theta >= thetaA);
            i2 = i2(1);
            
            dXdtheta = (obj.X(i2) - obj.X(i1))/(obj.theta(i2)-obj.theta(i1));
            dYdtheta = (obj.Y(i2) - obj.Y(i1))/(obj.theta(i2)-obj.theta(i1));
            
            XA = obj.X(i1) + dXdtheta*(thetaA - obj.theta(i1));
            YA = obj.Y(i1) + dYdtheta*(thetaA - obj.theta(i1));
            
            Phi = atan(dYdtheta/dXdtheta);
            
            ec = sin(Phi)*(X-XA) - cos(Phi)*(Y-YA);
            el = -cos(Phi)*(X-XA) - sin(Phi)*(Y-YA); 
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
