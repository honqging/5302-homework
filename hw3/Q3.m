x = 1:1:8;
y = 1:1:8;

% z(1~8)
z = [0,0,0,0,0,0,0,0;
	0,0,1,1,1,1,0,0;
	0,1,1,0,0,1,1,0;
	0,1,0,1,1,0,1,0;
	0,1,1,1,1,1,1,0;
	0,1,0,1,1,0,1,0;
	0,0,1,1,1,1,0,0;
	0,1,0,0,0,0,1,0]

xx = 1:0.02:8;
yy = 1:0.02:8;

% for x = 1.0, 1.02, 1.04, . . . , 7.98, 8.0
yys = zeros(8,351);

for i = 1:1:8
	pp = natspline(x,z(i,:));
	yys(i,:) = ppval(pp, xx);
end

% for y = 1.0, 1.02, 1.04, . . . , 7.98, 8.0
yys2 = zeros(351,351);
for j = 1:1:351
	p2 = natspline(y,yys(:,j)');
	yys2(:,j) = ppval(p2, yy)';
end

[x,y]=meshgrid(xx,yy);

% color map
 colormap(cool);
 mesh(x,y,yys2);
 colorbar;

% gray map
imagesc(yys2);
colormap(gray);
mesh(x,y,yys2);
colorbar;

title('2D spline');
xlabel('X');
ylabel('Y');
zlabel('Z');




