function MARK = reseed2_activ(MARK, asp, L, ney_t, nmy, ELEM2NODE_t, RHO_nodes, FUG_nodes, H_activ_nodes, PHASES_t, NUMMELT_t, TIMEMELT_t)

% auxiliary marker grid - new marker-nodes are the centers of the elements
[auxNODES, ELEMS, ~] =    generate_quad_grid(asp, nmy, 4, L);

% existing markers in auxiliary elements
[~, ~, MARK2AUX]   =    uv_quads(L, [asp*nmy; nmy], MARK.NODES);

% PHASE and RHO in auxiliary elements
sz                  =    [size(ELEMS,2),1];
Mweights            =    accumarray(MARK2AUX',ones(size(MARK2AUX)),sz); % counts how many markers in each aux elemet
PHASEsum            =    accumarray(MARK2AUX',MARK.PHASES,sz);
RHOsum              =    accumarray(MARK2AUX',MARK.RHO,sz);
FUGsum               =    accumarray(MARK2AUX',MARK.FUG,sz);
H_activsum         =    accumarray(MARK2AUX',MARK.H_activ,sz);
NUMMELTsum          =    accumarray(MARK2AUX',MARK.NUMMELT,sz);
TIMEMELTsum         =    accumarray(MARK2AUX',MARK.TIMEMELT,sz);

clear MARK MARK2AUX

MARK.PHASES         =    round(PHASEsum./Mweights)';
MARK.RHO            =    round(RHOsum./Mweights)';
MARK.FUG             =    round(FUGsum./Mweights)';
MARK.H_activ        =    round(H_activsum./Mweights)';
MARK.NUMMELT        =    round(NUMMELTsum./Mweights)';
MARK.TIMEMELT       =    round(TIMEMELTsum./Mweights)';
MARK.NODES          =    [mean(reshape(auxNODES(1,ELEMS),size(ELEMS,1),size(ELEMS,2)),1); mean(reshape(auxNODES(2,ELEMS),size(ELEMS,1),size(ELEMS,2)),1)];

% Fix NaNs (empty aux elements) by interpolating from the thermal grid
nanIND              =    isnan(MARK.RHO);
if sum(nanIND)>0
    nnodel_t             =    size(ELEM2NODE_t,1);
    [U, V, MARK2ELEM]    =    uv_quads(L, [asp*ney_t; ney_t], MARK.NODES(:,nanIND));
    SHP                  =    shp_quad([U; V], nnodel_t);
    MARK.RHO(nanIND)     =    sum(SHP.*RHO_nodes(ELEM2NODE_t(:,MARK2ELEM)));
    MARK.FUG(nanIND)      =    sum(SHP.*FUG_nodes(ELEM2NODE_t(:,MARK2ELEM)));
    MARK.H_activ(nanIND) =    sum(SHP.*H_activ_nodes(ELEM2NODE_t(:,MARK2ELEM)));
    MARK.PHASES(nanIND)  =    PHASES_t(MARK2ELEM);
    MARK.NUMMELT(nanIND) =    NUMMELT_t(MARK2ELEM);
    MARK.TIMEMELT(nanIND)=    TIMEMELT_t(MARK2ELEM);
end











