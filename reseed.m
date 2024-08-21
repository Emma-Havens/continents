function MARK = reseed(asp, L, ney_t, nmy, ELEM2NODE_t, RHO, RHO_nodes, PHASES_t)

dm                  =    2*L(2)/nmy;
Mx                  =    linspace(-(L(1)-dm), L(1)-dm, asp*nmy);
My                  =    linspace(-(L(2)-dm), L(2)-dm, nmy);
[MX, MY]            =    ndgrid(Mx,My);
MARK.NODES          =    [MX(:)'; MY(:)'];

nnodel_t            =    size(ELEM2NODE_t,1);
[U, V, MARK2ELEM]   =    uv_quads(L, [asp*ney_t; ney_t], MARK.NODES);
SHP                 =    shp_quad([U; V], nnodel_t);

MARK.RHO            =    sum(SHP.*RHO_nodes(ELEM2NODE_t(:,MARK2ELEM)));
MARK.RHO(MARK.RHO<mean(RHO))  = min(RHO);
MARK.RHO(MARK.RHO>=mean(RHO))  = max(RHO);

MARK.PHASES         =    PHASES_t(MARK2ELEM);