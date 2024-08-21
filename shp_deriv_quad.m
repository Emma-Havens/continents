function [N, dNdu] = shp_deriv_quad(ipx,nnodel)

nip               = size(ipx,1);
N                 = cell(nip,1);
dNdu              = cell(nip,1);

for i=1:nip
    eta1          = ipx(i,1);
    eta2          = ipx(i,2);
    switch nnodel
        case 4
            SHP   = 0.25*[(1-eta1)*(1-eta2);...
                (1+eta1)*(1-eta2);
                (1+eta1)*(1+eta2);
                (1-eta1)*(1+eta2)];
            
            DERIV = 0.25*[(-1+eta2) (1-eta2) (1+eta2) (-1-eta2);...  % w.r.t eta1
                (-1+eta1) (-1-eta1) (1+eta1) (1-eta1)];    % w.r.t eta2
            
        case 9
            N1X   = 0.5*eta1*(eta1-1);
            N2X   = -(eta1+1)*(eta1-1);
            N3X   = 0.5*eta1*(eta1+1);
            N1Y   = 0.5*eta2*(eta2-1);
            N2Y   = -(eta2+1)*(eta2-1);
            N3Y   = 0.5*eta2*(eta2+1);
            
            SHP   = [N1X*N1Y;...% 1 1
                N3X*N1Y;   % 3 1;
                N3X*N3Y;   % 3 3;
                N1X*N3Y;   % 1 3;
                N2X*N1Y;   % 2 1;
                N3X*N2Y;   % 3 2;
                N2X*N3Y;   % 2 3;
                N1X*N2Y;   % 1 2;
                N2X*N2Y];  % 2 2];
            
            DERIV = [(eta1-0.5)*N1Y (eta1+0.5)*N1Y (eta1+0.5)*N3Y...
                (eta1-0.5)*N3Y -2*eta1*N1Y (eta1+0.5)*N2Y...
                -2*eta1*N3Y (eta1-0.5)*N2Y -2*eta1*N2Y;...  % w.r.t eta1
                N1X*(eta2-0.5) N3X*(eta2-0.5) N3X*(eta2+0.5)...
                N1X*(eta2+0.5) N2X*(eta2-0.5) N3X*(-2)*eta2...
                N2X*(eta2+0.5) N1X*(-2)*eta2  N2X*(-2)*eta2];    % w.r.t eta2  %11 31 33 13 21 32 23 12 22;...
            
    end
    
    
    N{i}          = SHP;
    dNdu{i}       = DERIV';
    
end