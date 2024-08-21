function MARK = reseed2(MARK, asp, L, ney_t, nmy, ELEM2NODE_t, RHO, RHO_nodes, PHASES_t)

% auxiliary marker grid - new marker-nodes are the centers of the elements
[auxNODES, ELEMS, ~] =    generate_quad_grid(asp, nmy, 4, L);

% existing markers in auxiliary elements
[~, ~, MARK2AUX]   =    uv_quads(L, [asp*nmy; nmy], MARK.NODES);

% PHASE and RHO in auxiliary elements
sz                  =    [size(ELEMS,2),1];
Mweights            =    accumarray(MARK2AUX',ones(size(MARK2AUX)),sz); % counts how many markers in each aux elemet
PHASEsum            =    accumarray(MARK2AUX',MARK.PHASES,sz);
RHOsum              =    accumarray(MARK2AUX',MARK.RHO,sz);

temp_nummelt        =    MARK.NUMMELT;

clear MARK MARK2AUX

MARK.PHASES         =    round(PHASEsum./Mweights);
MARK.RHO            =    round(RHOsum./Mweights);
MARK.NODES          =    [mean(reshape(auxNODES(1,ELEMS),size(ELEMS,1),size(ELEMS,2)),1); mean(reshape(auxNODES(2,ELEMS),size(ELEMS,1),size(ELEMS,2)),1)];
MARK.NUMMELT        =    temp_nummelt;

% Fix NaNs (empty aux elements) by interpolating from the thermal grid
nanIND              =    isnan(MARK.RHO);
if sum(nanIND)>0
    nnodel_t            =    size(ELEM2NODE_t,1);
    [U, V, MARK2ELEM]   =    uv_quads(L, [asp*ney_t; ney_t], MARK.NODES(:,nanIND));
    SHP                 =    shp_quad([U; V], nnodel_t);
    MARK.RHO(nanIND)    =    sum(SHP.*RHO_nodes(ELEM2NODE_t(:,MARK2ELEM)));
    MARK.PHASES(nanIND) =    PHASES_t(MARK2ELEM);
end

MARK.RHO(MARK.RHO<mean(RHO))   = min(RHO);
MARK.RHO(MARK.RHO>=mean(RHO))  = max(RHO);

% Phases can only have values 1 or 2
MARK.PHASES(MARK.PHASES<1.5)   = 1;
MARK.PHASES(MARK.PHASES>=1.5)  = 2;










