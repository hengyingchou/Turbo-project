function [m,r,z] = input(p1,p2,p3,t,dt)

	%The purpose of this code is to find Bezier curve and find its length
 
	sum = 0;
 	X = 0:dt:t;
 	i=1;

	while(i ~= length(X))

		if length(X) ==1
    			break;
    		end

		drdt = (2*X(i)-2)*p1(2) + (2-4*X(i))*p2(2) + 2*X(i)*p3(2);
		dzdt = (2*X(i)-2)*p1(1) + (2-4*X(i))*p2(1) + 2*X(i)*p3(1);
		drdz = drdt/dzdt;

		dm = sqrt(drdz^2+1)*dzdt;
 
		%Using trapezoid rule to find the length of camber line 
		
		if i == 1 || i == length(X)
			sum = sum + dm;
		else
			sum = sum + 2*dm;    
		end

		i = i+1;

		end

			sum = sum *dt /2;

		m = sum;
		z = p1(1)*(1-t)^2 +p2(1)*2*t*(1-t)+p3(1)*t^2;
		r = p1(2)*(1-t)^2 +p2(2)*2*t*(1-t)+p3(2)*t^2;

	end
