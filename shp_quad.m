function SHP = shp_quad(uv,nnodel)

eta1         = uv(1,:);
eta2         = uv(2,:);

switch nnodel;
    case 4
        SHP = 0.25*[(1-eta1).*(1-eta2);...
            (1+eta1).*(1-eta2);
            (1+eta1).*(1+eta2);
            (1-eta1).*(1+eta2)];
        
    case 9
        N1X = 0.5*eta1.*(eta1-1);
        N2X = -(eta1+1).*(eta1-1);
        N3X = 0.5*eta1.*(eta1+1);
        N1Y = 0.5*eta2.*(eta2-1);
        N2Y = -(eta2+1).*(eta2-1);
        N3Y = 0.5*eta2.*(eta2+1);
        
        SHP = [N1X.*N1Y;
            N3X.*N1Y;
            N3X.*N3Y;
            N1X.*N3Y;
            N2X.*N1Y;
            N3X.*N2Y;
            N2X.*N3Y;
            N1X.*N2Y;
            N2X.*N2Y];
        
end