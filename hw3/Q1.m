knots = [-2 -1 0 1 2];
values=[1 1 0 1 1]; % for points y.
tt = -2:0.1:2;
p = polyfit(knots,values,4); % coefficients of interpolating polynomial
p
yyp = polyval(p,tt); % to evaluate the polynomial at all the points tt

pp = natspline(knots,values);% to compute the natural cubic spline
pp.coefs % to print out the coefficients of the spline
yys = ppval(pp,tt); % to evaluate the spline at many intermediate points 

plot(tt,yyp,'b-',tt,yys,'r--',knots,values,'ko');
% plot the polynomial interpolant (blue),
% the spline (red dashed), and the knots (black circles)
legend('interpolant','spline','knots','Location','Best')
% the last 2 arguments optimize the placement of the legend.
