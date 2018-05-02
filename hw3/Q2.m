k = 10;
knots = -1:1/k:1; % for points y.
values = 1./(1+100*knots.^2);
tt = -1:0.01:1;
ftValues = 1./(1+100*tt.^2); % all values 
p = polyfit(knots,values,4); % coefficients of interpolating polynomial
yyp = polyval(p,tt); % to evaluate the polynomial at all the points tt

pp = natspline(knots,values);% to compute the natural cubic spline
% pp.coefs % to print out the coefficients of the spline
yys = ppval(pp,tt); % to evaluate the spline at many intermediate points 

yyp(2);
plot(tt,ftValues,'g.',tt,yyp,'b-',tt,yys,'r--',knots,values,'ko');


maxError = 0;
for index = 1:1:201
	maxError = max(abs(ftValues(index) - yyp(index)), abs(ftValues(index) - yys(index)));
end

maxError;