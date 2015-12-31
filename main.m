% in this new version, phase shift is corrected by ray direction
clear
close all
%%
%parameters
lambda1=488; %nm
lambda2=592;
NA=1.4;
rMax=592;%lambda

zMin=-1184;
zMirror=0;
theta=0;
n=1.5;
rSize=50;
zSize=50;


Mag=100;
rMin=0;

theta1=0.01;
n1=1.5;

n2=0.13806+2.7024*1i;
n3=0.15078+3.5054*1i;





stage=1
%%
parameters.lambda1=lambda1; %nm
parameters.lambda2=lambda2;
parameters.NA=NA;
parameters.rMax=rMax;%lambda
parameters.rMin=rMin;
parameters.zMin=zMin;
parameters.zMirror=zMirror;
parameters.n=n;
parameters.n1=n1;
parameters.n2=n2;
parameters.n3=n3;
parameters.theta=theta;
parameters.rSize=rSize;
parameters.zSize=zSize;
theta2=asin(sin(theta1)*n1/n2);
theta22=asin(sin(theta1)*n1/n3);
parameters.phaseShift1=sin(theta2-theta1)/sin(theta2+theta1);
parameters.phaseShift2=sin(theta22-theta1)/sin(theta22+theta1);

stage=1
%%
    fieldEx=CalExcitationField(parameters,0);
    fieldEx2=CalExcitationField(parameters,1);
     fieldDep=CalDepletionField(parameters,0);
     fieldDep2=CalDepletionField(parameters,1);
%    
rCord=(rMax)/(rSize).*(-rSize:rSize);
zCord=zMin+(zMirror-zMin)/(zSize).*(0:2*zSize);
rCordNew=(rMax)/(rSize).*(-rSize:0.1:rSize);
zCordNew=zMin+(zMirror-zMin)/(zSize).*(0:0.1:2*zSize);
img1=abs(fieldEx').^2;
img2=abs(fieldEx2').^2;
img3=abs(fieldEx.'+fieldEx2.').^2;
 img4=abs(fieldDep').^2;
 img5=abs(fieldDep2').^2;
 img6=abs(fieldDep.'+fieldDep2.').^2;
[z0 r0]=meshgrid(zCord,rCord);
[zz rr]=meshgrid(zCordNew,rCordNew);
img1New=interp2(z0,r0,img1,zz,rr,'spline');
img4New=interp2(z0,r0,img4,zz,rr,'spline');
img3New=interp2(z0,r0,img3,zz,rr,'spline');
img6New=interp2(z0,r0,img6,zz,rr,'spline');
zMirrorLambda1=zMirror/lambda1;
mask=(zz+rr.*tan(theta)-zMirror<=0);


figure(1)
imagesc(zCordNew,rCordNew,img1New);axis image;
figure(2)
imagesc(zCord,rCord,abs(fieldEx2').^2);axis image;
figure(3)
imagesc(zCord,rCord,abs(fieldEx'+fieldEx2').^2);axis image;
figure(4)
imagesc(zCordNew,rCordNew,img3New.*mask)
xlabel('z/nm');ylabel('r/nm');axis image;
figure(5)
imagesc(zCordNew,rCordNew,img6New.*mask)
xlabel('z/nm');ylabel('r/nm');axis image;
