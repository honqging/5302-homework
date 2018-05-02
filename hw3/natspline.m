function pp = natspline(x,y,conds)
%CSAPE Cubic spline interpolation with various end conditions.
%   Default changed to be the 'natural' cubic spline (2nd deriv's == 0 at ends)
%
%   PP  = CSAPE(X,Y)  returns the cubic spline interpolant (in ppform) to the 
%   given data (X,Y) using Lagrange end conditions (the default in table below).
%   The interpolant matches, at the data site X(j), the given data value
%   Y(:,j), j=1:length(X). The data values may be scalars, vectors.
%
%   PP  = CSAPE(X,Y,CONDS)  uses the end conditions specified by CONDS, with
%   corresponding end condition values  endcondvals .
%   If there are two more data values than data sites, then the first (last) 
%   data value is taken as the value for the left (right) end condition, i.e.,
%   endcondvals = Y(:,[1 end]).
%   Otherwise, default values are used.
%
%   CONDS may be a *string* whose first character matches one of the
%   following: 'complete' or 'clamped', 'periodic',
%   'second', 'variational', with the following meanings:
%
%   'complete'    : match endslopes to the slope of the cubic that
%                   matches the first four data at the respective end.
%   'not-a-knot'  : no longer supported by this function.
%   'periodic'    : match first and second derivatives at first data
%                   point with those at last data point
%                   (ignoring given end condition values if any)
%   'second'      : match end second derivatives (as given,
%                   with default [0 0], i.e., as in variational)
%   'variational' : set end second derivatives equal to zero
%                   (ignoring given end condition values if any)
%   The *default* : natural cubic spline (like 'second' w/ zero end conditions)
%
%   By giving CONDS as a 1-by-2 matrix instead, it is possible to
%   specify *different* conditions at the two endpoints, namely
%   CONDS(i) with value endcondvals(:,i), with i=1 (i=2) referring to the
%   left (right) endpoint.
%
%   CONDS(i)=j  means that the j-th derivative is being specified to
%   be endcondvals(:,i) , j=1,2.  CONDS(1)=0=CONDS(2)  means periodic end
%   conditions.
%
%   If CONDS(i) is not specified or is different from 0, 1 or 2, then
%   the default value for CONDS(i) is  1  and the default value of
%   endcondvals(:,i) is taken.  If no end condition values are specified,
%   then the default value for endcondvals(:,i) is taken to be
%
%    deriv. of cubic interpolant to nearest four points, if   CONDS(i)=1;
%                     0                                  if   CONDS(i)=2.
%
%   For example,
%
%      x = linspace(0,2*pi,9);    %% sample every 45 degrees.
%      pp = natspline( x, [1 sin(x) 0], [1 2] );
%      xx=linspace(0,2*pi,50);    %% many samples for drawing curve.
%      plot(xx,ppval(pp,xx));
%
%   gives a good approximation to the sine function on the interval [0 .. 2*pi]
%   (matching its slope 1 at the left endpoint, x(1) = 0, and its second 
%   derivative 0 at the right endpoint, x(9) = 2*pi, in addition to its value
%   at every x(i), i=1:9).
%
%   The following plots a circle:
%
%      x = [0:.1:4]; pp=natspline( [0:4], [1 0 -1 0 1;0 1 0 -1 0], 'periodic');
%      y = ppval(pp,x);  plot(y(1,:),y(2,:));  axis equal;
%
%   Copyright 1987-2003 C. de Boor and The MathWorks, Inc.
%   $Revision: 1.23 $   editted for use in UofM class
%     Generate the cubic spline interpolant in ppform.

if nargin<3, conds = [2 2]; end

%     Generate the cubic spline interpolant in ppform.
% The fourth argument still permitted here for backward compatibility

% [xi,yi,sizeval,endvals] = chckxywp(x,y,0);
xi=x(:);
sizeval=length(x);
if sum(size(x))-1 ~= sizeval;
   error('natspline: input data must be simple vectors.');
end;
if length(y) == sizeval;
   yi=y(:);
   endvals=[];
elseif length(y)==sizeval+2;
   yi = y(2:sizeval+1);
   yi=yi(:);
   endvals=[y(1);y(sizeval+2)];
else
   error('natspline: mismatch in vector lengths in input arguments');
end; 

if ~isempty(endvals), valconds = endvals.'; end
[yn,yd] = size(yi); dd = ones(1,yd);
dx = diff(xi); divdif = diff(yi)./dx(:,dd);

[n,yd] = size(yi); dd = ones(1,yd);

valsnotgiven=0;
if ~exist('valconds','var'), valsnotgiven=1;  valconds = zeros(yd,2); end
if ischar(conds)
   if     conds(1)=='c', conds = [1 1];
   %elseif conds(1)=='n', conds = [1 1]; 
   elseif conds(1)=='p', conds = [0 0];
   elseif conds(1)=='s', conds = [2 2];
   elseif conds(1)=='v', conds = [2 2]; valconds = zeros(yd,2);
   else, error('SPLINES:CSAPE:unknownends',...
              ['Unknown end condition *',conds,'* specified.'])
   end
end

   % set up the linear system for solving for the slopes at XI.
dx = diff(xi); divdif = diff(yi)./dx(:,dd);
c = spdiags([ [dx(2:n-1,1);0;0] ...
            2*[0;dx(2:n-1,1)+dx(1:n-2,1);0] ...
              [0;0;dx(1:n-2,1)] ], [-1 0 1], n, n);
b = zeros(n,yd);
b(2:n-1,:)=3*(dx(2:n-1,dd).*divdif(1:n-2,:)+dx(1:n-2,dd).*divdif(2:n-1,:));
if ~any(conds)
   c(1,1)=1; c(1,n)=-1;
elseif conds(1)==2
   c(1,1:2)=[2 1]; b(1,:)=3*divdif(1,:)-(dx(1)/2)*valconds(:,1).';
else
   c(1,1:2) = [1 0]; b(1,:) = valconds(:,1).';
   if ((valsnotgiven)||(conds(1)~=1))  % if endslope was not supplied,
                                   % get it by local interpolation
     b(1,:)=divdif(1,:);
     if n>2, ddf=(divdif(2,:)-divdif(1,:))/(xi(3)-xi(1));
       b(1,:) = b(1,:)-ddf*dx(1); end
     if n>3, ddf2=(divdif(3,:)-divdif(2,:))/(xi(4)-xi(2));
       b(1,:)=b(1,:)+(ddf2-ddf)*(dx(1)*(xi(3)-xi(1)))/(xi(4)-xi(1)); end
   end
end
if ~any(conds)
   c(n,1:2)=dx(n-1)*[2 1]; c(n,n-1:n)= c(n,n-1:n)+dx(1)*[1 2];
   b(n,:) = 3*(dx(n-1)*divdif(1,:) + dx(1)*divdif(n-1,:));
elseif conds(2)==2
   c(n,n-1:n)=[1 2]; b(n,:)=3*divdif(n-1,:)+(dx(n-1)/2)*valconds(:,2).';
else
   c(n,n-1:n) = [0 1]; b(n,:) = valconds(:,2).';
   if ((valsnotgiven)||(conds(2)~=1))  % if endslope was not supplied,
                                   % get it by local interpolation
      b(n,:)=divdif(n-1,:);
      if n>2, ddf=(divdif(n-1,:)-divdif(n-2,:))/(xi(n)-xi(n-2));
        b(n,:) = b(n,:)+ddf*dx(n-1); end
      if n>3, ddf2=(divdif(n-2,:)-divdif(n-3,:))/(xi(n-1)-xi(n-3));
        b(n,:)=b(n,:)+(ddf-ddf2)*(dx(n-1)*(xi(n)-xi(n-2)))/(xi(n)-xi(n-3));
      end
   end
end

  % solve for the slopes ..  (protect current spparms setting)
mmdflag = spparms('autommd');
spparms('autommd',0); % suppress pivoting
s=c\b;
spparms('autommd',mmdflag);

  %                          .. and convert to ppform
c4 = (s(1:n-1,:)+s(2:n,:)-2*divdif(1:n-1,:))./dx(:,dd);
c3 = (divdif(1:n-1,:)-s(1:n-1,:))./dx(:,dd) - c4;
pp =  mkpp(xi.', ...
   reshape([(c4./dx(:,dd)).' c3.' s(1:n-1,:).' yi(1:n-1,:).'],(n-1)*yd,4),yd);
if length(sizeval)>1, pp = fnchg(pp,'dz',sizeval); end

if ~isfield(pp,'coefs');pp.coefs=pp.P;end;
