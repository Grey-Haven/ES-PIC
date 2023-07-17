function F_p = gather_field_2D(F,loc,x,y)
    
    index_offset = 1;
    a_x = x(1);
    a_y = y(1);
    
    dx = x(2)-x(1);
    dy = y(2)-y(1);

    lc_x = index_offset+(loc(:,1)-a_x)/dx;
    lc_y = index_offset+(loc(:,2)-a_y)/dy;

    is = floor(lc_x);
    js = floor(lc_y);

    x_node = x(is)';
    y_node = y(js)';

    d1 = (loc(:,1) - x_node)./dx;
    d2 = (loc(:,2) - y_node)./dy;

    idx1 = sub2ind(size(F),js,is);
    idx2 = sub2ind(size(F),js+1,is);
    idx3 = sub2ind(size(F),js,is+1);
    idx4 = sub2ind(size(F),js+1,is+1);

    F1 = F(idx1).*(1-d1).*(1-d2);
    F2 = F(idx2).*d1.*(1-d2);
    F3 = F(idx3).*(1-d1).*d2;
    F4 = F(idx4).*d1.*d2;
    
    F_p = F1+F2+F3+F4;
end