function rho = scatter_charge_2D(x_locs, y_locs, qs, ws, x, y, Nx, Ny)

    dx = x(2)-x(1);
    dy = y(2)-y(1);

    index_offset = 1;

    a_x = x(1);
    a_y = y(1);
    
    lc_x = index_offset+(x_locs(:) - a_x)./dx;
    lc_y = index_offset+(y_locs(:) - a_y)./dy;
    is = floor(lc_x);
    js = floor(lc_y);
    wzs = qs.*ws;

    x_nodes = x(is)';
    y_nodes = y(js)';

    fxs = (x_locs(:) - x_nodes)/dx;
    fys = (y_locs(:) - y_nodes)/dy;

    assert(all(0 <= fxs & fxs <= 1))

    rho5 = accumarray([js,is],(1-fxs).*(1-fys).*wzs,[Nx,Ny]);
    rho6 = accumarray([js,is+1],(1-fxs).*fys.*wzs,[Nx,Ny]);
    rho7 = accumarray([js+1,is],fxs.*(1-fys).*wzs,[Nx,Ny]);
    rho8 = accumarray([js+1,is+1],fxs.*fys.*wzs,[Nx,Ny]);

    rho = (rho5(:,:)+rho6(:,:)+rho7(:,:)+rho8(:,:)) ./ (dx*dy);
end