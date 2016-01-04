% in this new version, phase shift is corrected by ray direction
clear
close all
%%
%parameters
lambda1=488; %extation wavelength, in the unit of nm
lambda2=592;%depletion wavelength, in the unit of nm
NA=1.4;% numerical aperture of objective
rMax=592;%simulation area, in unit of nm

zMin=-1184;%simulation area, in unit of nm
zMirror=0;%the distance from focal point to mirror, in unit of nm
theta=0;%the angle mirror tilt
n=1.5;%refractive index of image plane
rSize=50;% grids of simulation
zSize=50;%grids of simulation

Mag=100;%magnification. only used in amplitude
rMin=0;%not used

theta1=0.01;%not used, for test only
n1=1.5;%not used, for test only

n2=0.13806+2.7024*1i;%refractive index of silver, for excitation 
n3=0.15078+3.5054*1i;%refractive index of silver, for depletion 





stage=1
%%
%group the parameters to be used in functuions
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

stage=2
%%
    fieldEx=CalExcitationField(parameters,0);%calculate the incident excitation field
    fieldEx2=CalExcitationField(parameters,1);%calculate the reflected excitation field
     fieldDep=CalDepletionField(parameters,0);%calculate the incident depletion field
     fieldDep2=CalDepletionField(parameters,1);%calculate the reflected depletion field
%    
rCord=(rMax)/(rSize).*(-rSize:rSize);%coordinatin along r
zCord=zMin+(zMirror-zMin)/(zSize).*(0:2*zSize);%coordinatin along z
rCordNew=(rMax)/(rSize).*(-rSize:0.1:rSize);%coordinatin along r after interpletion
zCordNew=zMin+(zMirror-zMin)/(zSize).*(0:0.1:2*zSize);%coordinatin along z after interpletion
%calculate the intensity distributions
img1=abs(fieldEx').^2;%incident excitation
img2=abs(fieldEx2').^2;%reflected excitation
img3=abs(fieldEx.'+fieldEx2.').^2;%total intensity
 img4=abs(fieldDep').^2;%incident depletion
 img5=abs(fieldDep2').^2;%reflected depletion
 img6=abs(fieldDep.'+fieldDep2.').^2;%overall depletion
%interplate to get  higher resolution images
 [z0 r0]=meshgrid(zCord,rCord);%old coordination
[zz rr]=meshgrid(zCordNew,rCordNew);%new coordination
img1New=interp2(z0,r0,img1,zz,rr,'spline');
img4New=interp2(z0,r0,img4,zz,rr,'spline');
img3New=interp2(z0,r0,img3,zz,rr,'spline');
img6New=interp2(z0,r0,img6,zz,rr,'spline');
zMirrorLambda1=zMirror/lambda1;
mask=(zz+rr.*tan(theta)-zMirror<=0);%use a mask to hide the area behind the mirror
%%
%draw the pictures
% figure(1)
% imagesc(zCordNew,rCordNew,img1New);axis image;
% figure(2)
% imagesc(zCord,rCord,abs(fieldEx2').^2);axis image;
% figure(3)
% imagesc(zCord,rCord,abs(fieldEx'+fieldEx2').^2);axis image;
figure(4)
imagesc(zCordNew,rCordNew,img3New.*mask)
xlabel('z/nm');ylabel('r/nm');axis image;
figure(5)
imagesc(zCordNew,rCordNew,img6New.*mask)
xlabel('z/nm');ylabel('r/nm');axis image;
