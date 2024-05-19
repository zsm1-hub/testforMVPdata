function [xddd2,addd2]=avg_alongtrack(lon,lat,xxx,a)
dist1=spheric_dist(lat(2:end),lat(1:end-1),lon(2:end),lon(1:end-1));
xdd1(1)=0;
for len=2:length(lat)
    xdd1(len)=xdd1(1)+sum(dist1(1:len-1));
end

for i=1:length(xxx)-1
    dot1=find(xdd1>=xxx(i));
    xddd1=xdd1(dot1);
    addd1=a(dot1);
    dot2=find(xddd1<=xxx(i+1));
    xddd2=xddd1(dot2);
    addd2(i)=nanmean(addd1(dot2));
end
return