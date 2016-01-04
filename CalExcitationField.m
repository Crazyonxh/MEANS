function filedEx=CalExcitationField(parameters,reflect)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

lambda1=parameters.lambda1; %nm
NA=parameters.NA;
rMax=parameters.rMax;%lambda
rMin=-parameters.rMax;
zMin=parameters.zMin;
zMirror=parameters.zMirror;
n=parameters.n;
rSize=parameters.rSize;
zSize=parameters.zSize;
thetaMirror=parameters.theta;
n1=parameters.n1;
n2=parameters.n2;
rMaxLambda1=rMax/lambda1;%lambda
rMinLambda1=rMin/lambda1;
zMinLambda1=zMin/lambda1;
zMirrorLambda1=zMirror/lambda1;
rCord=(rMaxLambda1)/(rSize).*(-rSize:rSize);
zCord=zMinLambda1+(zMirrorLambda1-zMinLambda1)/(zSize).*(0:2*zSize);
%%
k = 2*pi*n/lambda1;
alpha = asin(NA/n);
filedEx=zeros(2*zSize+1,rSize+1);
u=zeros(2*zSize+1,2*rSize+1);
v=zeros(2*zSize+1,2*rSize+1);


%%
if reflect==0
    for i=(1:(2*zSize+1))
        for j=1:(2*rSize+1)
          u(i,j)=4*k*zCord(i)*lambda1*(sin(alpha/2)^2);
          v(i,j)=k*rCord(j)*lambda1*sin(alpha); 
        end;
    end;
    for i=1:(2*zSize+1)
        for j=1:(2*rSize+1)
        Koi = -2*pi*1i/(lambda1)*exp(1i*u(i,j)/(4*(sin(alpha/2)^2)));
        intgrand = @(theta) (sqrt(cos(theta))) .*  (sin(theta)) .*   (exp((-1i*u(i,j)/2)* (sin(theta/2).^2) / (sin(alpha/2)^2)))  .*  (besselj(0, sin(theta)/sin(alpha).*v(i,j)));%Min Gu P149
        I0 = integral(@(theta)intgrand (theta),0,alpha);  
      filedEx(i,j)=Koi.*I0;
        end;
    end;
end
if reflect==1
    for i=(1:(2*zSize+1))
        for j=1:(2*rSize+1)
            zz=(tan(thetaMirror).^2-1)/(1+tan(thetaMirror).^2).*zCord(i)+(-tan(thetaMirror)*2)/(1+tan(thetaMirror).^2).*rCord(j)+(2)/(1+tan(thetaMirror).^2).*zMirrorLambda1;
            rr=(-tan(thetaMirror)*2)/(1+tan(thetaMirror).^2).*zCord(i)-(tan(thetaMirror).^2-1)/(1+tan(thetaMirror).^2).*rCord(j)+(2*tan(thetaMirror))/(1+tan(thetaMirror).^2).*zMirrorLambda1;
         u(i,j)=4*k*zz*lambda1*(sin(alpha/2)^2);
         v(i,j)=k*rr*lambda1*sin(alpha);    
        end;
    end;
    for i=1:(2*zSize+1)
        for j=1:(2*rSize+1)
        Koi = -2*pi*1i/(lambda1)*exp(1i*u(i,j)/(4*(sin(alpha/2)^2)));
       intgrand = @(theta) (sqrt(cos(theta))) .*  (sin(theta)) .*   (exp((-1i*u(i,j)/2)* (sin(theta/2).^2) / (sin(alpha/2)^2)))  .*  (besselj(0, sin(theta)/sin(alpha).*v(i,j))).* phaseShiftCalculation( n1,n2,theta);%Min Gu P149
       %  intgrand = @(theta) (sqrt(cos(theta))) .*  (sin(theta)) .*   (exp((-1i*u(i,j)/2)* (sin(theta/2).^2) / (sin(alpha/2)^2)))  .*  (besselj(0, sin(theta)/sin(alpha).*v(i,j)));%Min Gu P149
        I0 = integral(@(theta)intgrand (theta),0,alpha);  
      filedEx(i,j)=Koi.*I0;
        end;
    end;
end;


end
