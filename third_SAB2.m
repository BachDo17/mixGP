function [Pf,PDF]=third_SAB2(muy1,muy2,muy3,y)
% INPUTS
% muy1,muy2,muy3 = first three moment of the limit state function eveluated
% at the mean vector of input variables
% y = threshold of the limit state function

% OUTPUT
% Pf: the failure probability
% PDF: of the limit state function
if muy3 ==0
  a = 0;
  b = 0;
else
  a = 2*muy2^3/muy3^2;
  b = 0.5*muy3/muy2;
end

% Saddlepoint
Sadd = (y-muy1)/(-b*muy1+muy2+y*b);
% Compute r and v
diff0 = (muy1-2*a*b)*Sadd+0.5*(muy2-2*a*b^2)*Sadd^2-a*log((1-b*Sadd)^2);
diff2 = muy2-2*a*b^2+(2*a*b^2)/(1-b*Sadd)^2;
r = sign(Sadd)*sqrt(2*(Sadd*y-diff0));
v= Sadd*sqrt(diff2);
% PDF
PDF = 1/(sqrt(2*pi*diff2))*exp(diff0-Sadd*y+0.5*Sadd^2*diff2)*exp(-0.5*Sadd^2*diff2);
% CDF
if y==muy1
    Pf = 0.5+muy3/(6*sqrt(2*pi)*(muy2)^1.5);
else
    Pf = normcdf(r) + normpdf(r)*(1/r - 1/v);
end

end