 The key element in this code is in CalExcitationField.m
       intgrand = @(theta) (sqrt(cos(theta))) .*  (sin(theta)) .*   (exp((-1i*u(i,j)/2)* (sin(theta/2).^2) / (sin(alpha/2)^2)))  .*  (besselj(0, sin(theta)/sin(alpha).*v(i,j))).* phaseShiftCalculation( n1,n2,theta);
       from Min Gu P 149, where a conjugation is used to be consistant with Wolf 
       The phase calculation is from Wolf
