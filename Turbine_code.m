clc;
clear;
close;

%Basic parameter for design

 inlet = 660;
 omega = 800/0.98; 
 r = 0.2;
 dr = 0.78/3;
 P1 = [ 0; 0.98]; P2 = [0; 0.2]; P3 = [1.18412 ; 2];P4 = [1.3682;2];

	for i = 1 : 1 : 4
    		%Bezir Curve

    		if i==1

    			p1 = P2;
    			p2(1) = P4(1);  
    			p2(2) = P2(2);
    			p3 = P4;

    		elseif i ==4 
	    		p1 =[0 ; 0.98];
    			p2 = [1.18412;0.98];
    			p3 = P3;
	    		else
    			p1(1) = 0;     
    			p1(2)= (0.78/3)*(i-1) + 0.2;
	    		p3(1) =1.3682 - (0.2301/3)*(i-1);
    			P3(2) = 2;
    			p2(1) = p3(1);
    			p2(2) = p1(2);
    		end

	e = 0 : 0.001 : 0.9999999;
	for k = 1 : length(e)
         	[a,b,c]  = input(p1,p2,p3,e(k),0.001);
          	m(k) = a;
          	R(k) = b;
          	z(k) = c;
    	end

	l = max(m);
    	cam(i) = l;

% Initial theta for calcualting mental angle 

   	r = 0.2+(i-1)*dr; 
   	rotational_speed = omega*r;
   	dthetadm = rotational_speed/(inlet*r);
   	rotational_speed2 = omega*2;  
   	dthetadm2 = (omega*2 - 660)/(2*660*cos(pi/3));
   	Matrix =[0 0 0 1; 0 0 1 0;l^3 l^2 l 1; 3*l^2 2*l 1 0];
   	C = [ 0 ; dthetadm; (50-(i-1)*5/3)*pi/180; 0];
   	B = inv(Matrix)*C;

   	for j = 1 : length(m)
        	W(j) = B(1)*m(j)^3 + B(2)*m(j)^2 + B(3)*m(j) + B(4);
        	beta(j) = atan((3*B(1)*m(j)^2 + 2*B(2)*m(j)+B(3))*R(j));
   	end

   	for j = 1 : length(m)
   		x(j) = R(j)*cos(W(j));
   		y(j)=  R(j)*sin(W(j));
   		dydx(j) = tan(W(j));
   	end

% Thickness distribution
 
	%calculate maximum thickness
	t = 0.02*l;
 
	%calcuate thickness distribution
 	for j = 1 :100

		p = m(j)/m(100);
		yt(j) = t *5*( (0.2969*sqrt(p) - 0.1260*p - 0.3516*p^2 + 0.2843*p^3 - 0.1015*p^4));
 	
	end

	count = 0;
		for j = 1:100

			if (yt(j) == max(yt))
				count = j;
			end
	end

	rate = t/max(yt);
	yt = yt *rate;

	for  j = count : length(m)-count

		yt(j) = yt(count);
	end

	for j = 1 : count

		Q(j) = yt(j);
	end

	for j = 1 : count
 
		EMU(j) = Q(count-j+1) ;
	
	end

	yt = [yt,EMU];


% Making Airfoil

	for j = 1 : length(m)   
 
   		xu(j) =  x(j)+yt(j) * -sin(dydx(j));
   		yu(j) = y(j) + yt(j) * cos(dydx(j));
   		xl(j) = x(j)- yt(j) * -sin(dydx(j));
   		yl(j) = y(j)-yt(j) * cos(dydx(j));

	end

	filename = sprintf('U%d.txt',i);
	filename2 = sprintf('L%d.txt',i);
    	fileID = fopen(filename,'w');    
    	fileID2 = fopen(filename2,'w');

    	for j = 1 : length(m)
		fprintf(fileID,'%f\t%f\t%f\r\n',xu(j)*304.7,yu(j)*304.7,z(j)*304.7);
		fprintf(fileID2,'%f\t%f\t%f\r\n',xl(j)*304.7,yl(j)*304.7,z(j)*304.7);
    	end

	fclose(fileID);
	fclose(fileID2);

% Making hub

	if i ==1

		filID1 = fopen('Hub.txt','w'); 
		
		p1 = [0;0.201];
		p2 = [1.367;0.201];
		p3 = [1.367;2];
		
		TT = 0: 0.01 : 1;
		for j = 1 : length(TT)
			Z(j) = (1-TT(j))^2*p1(1)+ 2*TT(j)*(1-TT(j))*p2(1)+TT(j)^2*p3(1);
			R2(j) = (1-TT(j))^2*p1(2)+ 2*TT(j)*(1-TT(j))*p2(2)+TT(j)^2*p3(2);
			fprintf(filID1,'%f\t%f\t%f\r\n',0,R2(j)*304.7,Z(j)*304.7);
		end

		fclose(filID1);
	end
% Making Duct

	if i == 4

    		ID = fopen('Duct.txt','w');

	for j = 1 : 1 : length(z)

    		if j ==1
    
			p1(1) = 0;
    			p2(1) = 1.0;
    
		elseif j== length(z)
     
			p1(j) = 1.1668;
     			p2(j) = 2.0;
    
		else
        
		a1 = 2*z(j) - 2*z(j-1);
	        b1 = 2*R(j) - 2*R(j-1);
        	c1 = z(j)^2 + R(j)^2 -z(j-1)^2 -R(j-1)^2 ;
	        a2 = 2*z(j+1) - 2*z(j);
        	b2 = 2*R(j+1) - 2*R(j);
 	       	c2 = z(j+1)^2 + R(j+1)^2 -z(j)^2 -R(j)^2 ;
       		
		M = [a1 b1;a2 b2];     
       		U = [c1;c2]; 
       		OO = inv(M)*U; 
       
		x = OO(1); 
       		y = OO(2);
       		l = sqrt((x-z(j))^2 + (y-R(i))^2);
       
		x_d = (x-z(j))/l;
       		y_d = (y-R(j))/l;
       
		p1(j) = 0.02*x_d + z(j);
       		p2(j) = 0.02*y_d + R(j);
		end

	fprintf(ID,'%f\t%f\t%f\r\n',0,p2(j)*304.7,p1(j)*304.7);
	
	end

	fclose(ID);	
	end

	xu = [];
	zu = [];
	xl = [];
	zl = [];	
	yt = [];
   	x = [];
   	y = [];
	m = [];
	r = [];
	z = [];
	W = [];

	end
