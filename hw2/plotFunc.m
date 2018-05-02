x1 = 59.4727;
x2 = 42.6599;
    a = x1;

    g = 9.8065;
    r = 0.3;
    v = 15;

   
x=linspace(0,11); 

y = x * (g/(r * v * cosd(a)) + tand(a)) + g/(r^2) * log(1 - r * x/(v * cosd(a)));


plot(x, y, 'r')
hold on

a = x2;
y2 = x * (g/(r * v * cosd(a)) + tand(a)) + g/(r^2) * log(1 - r * x/(v * cosd(a)));

plot(x, y2, 'b')
legend('degree: 59.4727', 'degree: 42.6599');