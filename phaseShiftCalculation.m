function PhaseShift = phaseShiftCalculation( n1,n2,theta )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
theta2=asin(sin(theta)*n1/n2);
PhaseShift=sin(theta2-theta)/sin(theta2+theta);
end

