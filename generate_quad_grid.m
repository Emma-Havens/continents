function [GCOORD, ELEM2NODE, ind] = generate_quad_grid(asp, ney, nnodel, L)

asp = round(asp);

switch nnodel
    case 9
        nex          =  asp*ney;
        nel          =  nex*ney;
        
        nx           =  2*nex+1;
        ny           =  2*ney+1;
        x            =  linspace(-L(1),L(1),nx);
        y            =  linspace(-L(2),L(2),ny);
        
        [X Y]        =  ndgrid(x, y);
        GCOORD       =  [X(:)';Y(:)'];
        
        ind          =  reshape(1:nx*ny, nx, ny);
        ELEM2NODE    =  zeros(nnodel,nel);
        ri           =  [0 2 2 0 1 2 1 0 1];
        rj           =  [0 0 2 2 0 1 2 1 1];
        for k=1:9
            tmp      =  ind(1+ri(k):2:end-2+ri(k),1+rj(k):2:end-2+rj(k)) ;
            ELEM2NODE(k,:) = tmp(:);
        end
        ELEM2NODE    =  int32(ELEM2NODE);
    case 4
        nex          =  asp*ney;
        nel          =  nex*ney;
        
        nx           =  nex+1;
        ny           =  ney+1;
        x            =  linspace(-L(1),L(1),nx);
        y            =  linspace(-L(2),L(2),ny);
        
        [X Y]        =  ndgrid(x, y);
        GCOORD       =  [X(:)';Y(:)'];
        
        ind          =  reshape(1:nx*ny,nx, ny);
        ELEM2NODE    =  zeros(nnodel,nel);
        ri           =  [0 1 1 0];
        rj           =  [0 0 1 1];
        for k=1:4
            tmp      =  ind(1+ri(k):1:end-1+ri(k),1+rj(k):1:end-1+rj(k)) ;
            ELEM2NODE(k,:) = tmp(:);
        end
        ELEM2NODE  	 =  int32(ELEM2NODE);
end