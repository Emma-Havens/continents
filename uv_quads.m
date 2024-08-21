function [U, V, IND] =  uv_quads(L,nel, GCOORD)

dx = 2*L./nel;

% SHIFTING THE DOMAIN TO LIE IN THE FIRST QUADRANT
GCOORD(1,:)   = GCOORD(1,:) + L(1);
GCOORD(2,:)   = GCOORD(2,:) + L(2);

% FINDING TO WHICH ELEMENTS POINTS GCOORD BELONG
indx          = ceil(GCOORD(1,:)./dx(1));
indx(indx<1)  = 1; indx(indx>nel(1)) = nel(1);
indy          = ceil(GCOORD(2,:)./dx(2));
indy(indy<1)  = 1; indy(indy>nel(2)) = nel(2);

IND           = indx + (indy-1).*nel(1);

% SCALING POSITIONS GCOORD TO LOCAL ELEMENT POSITIONS
U             = 2*(GCOORD(1,:) - (indx-1+1/2)*dx(1))./dx(1);
V             = 2*(GCOORD(2,:) - (indy-1+1/2)*dx(2))./dx(2);