yTarget = 4;
a = 0;





i = 0;
j = 45;
n = 0;

while(j - i > 0.01)
    n = n+1;
    m = i + (j - i)/2;
    if((water(m) - yTarget) * (water(j) - yTarget) < 0)
        i = m;
    else
        j = m;
    end
end

i;
j;
ns = water(i);
n2s = water(j);