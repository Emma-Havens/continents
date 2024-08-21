function [Vel, Pressure] = mechanical2d_simple(ELEM2NODE, PHASES, GCOORD, Mu, Rho, Ra, G, bc_ind, bc_val, nip)
                                             
div_max_uz     =    1e-10; div_max  = realmax;
uz_iter        =    0; uz_iter_max  =       5;

nel               = size(ELEM2NODE,2);
nnodel            = size(ELEM2NODE,1);
nnod              = size(GCOORD,2);
ndim              = size(GCOORD,1);

%==========================================================================
% CONSTANTS
%==========================================================================
nedof       = nnodel*ndim;
sdof        = 2*nnod;
np          = 3;

DEV   = [ 4/3 -2/3 0;...
         -2/3  4/3 0;...
            0    0 1];

PF = 1e3*max(Mu(:));

[IP_X, IP_w] = ip_line(nip);
[xx, yy] = meshgrid(IP_X,IP_X);
IP_X = [xx(:) yy(:)];
IP_w = IP_w'*IP_w;
IP_w = IP_w(:);
nip  = length(IP_w);
[N,     dNdu] = shp_deriv_quad(IP_X,nnodel);

%==========================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==========================================================================
A_all       = zeros(sdof,sdof);
Q_all       = zeros(np*nel,sdof);
invM_all    = zeros(np*nel,np*nel);
Rhs_all     = zeros(sdof,1);

A_elem      = zeros(nedof,nedof);
Q_elem      = zeros(nedof,np);
M_elem      = zeros(np,np);
Rhs_elem    = zeros(ndim,nnodel);

B           = zeros(nedof,ndim*(ndim+1)/2);
P           = ones(np);
Pb          = ones(np,1);

for iel = 1:nel
    %==============================================================
    % FETCH DATA OF ELEMENT
    %==============================================================
    ECOORD_X  = GCOORD(:,ELEM2NODE(:,iel));
    EMu       = Mu(PHASES(iel));
    ERho      = Rho(PHASES(iel));
    
    %==============================================================
    % INTEGRATION LOOP
    %==============================================================
    A_elem(:) = 0;
    Q_elem(:) = 0;
    M_elem(:) = 0;
    Rhs_elem(:) = 0;
    
    P(2:3,:) = ECOORD_X(:,1:3);
    for ip=1:nip
        %==========================================================
        % LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
        %==========================================================
        Ni          =       N{ip};
        dNdui       =    dNdu{ip};
        Pb(2:3)     = ECOORD_X*Ni;
        Pi          = P\Pb;
        
        %==========================================================
        % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
        %==========================================================
        J           = ECOORD_X*dNdui;
        detJ        = det(J);
        invJ        = inv(J);
        
        %==========================================================
        % DERIVATIVES wrt GLOBAL COORDINATES
        %==========================================================
        dNdX        = dNdui*invJ;
        
        %==========================================================
        % NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==========================================================
        weight       = IP_w(ip)*detJ;
        B(1:2:end,1) = dNdX(:,1);
        B(2:2:end,2) = dNdX(:,2);
        B(1:2:end,3) = dNdX(:,2);
        B(2:2:end,3) = dNdX(:,1);
        Bvol         = dNdX';
        
        A_elem       = A_elem + weight*EMu*(B*DEV*B');
        Q_elem       = Q_elem - weight*Bvol(:)*Pi';
        M_elem       = M_elem + weight*Pi*Pi';
        Rhs_elem(1,:)= Rhs_elem(1,:) + weight*ERho*G(1)*Ni';
        Rhs_elem(2,:)= Rhs_elem(2,:) + weight*ERho*G(2)*Ni';
    end
    %==============================================================
    % STATIC CONDENSATION
    %==============================================================
    invM_elem = inv(M_elem);
    A_elem    = A_elem + PF*Q_elem*invM_elem*Q_elem';
    
    %==============================================================
    % WRITE DATA INTO GLOBAL STORAGE
    %==============================================================
    ind_el           = zeros(1,nedof);
    ind_el(1:2:end-1)= 2*(ELEM2NODE(:,iel)-1)+1;
    ind_el(2:2:end)  = 2*(ELEM2NODE(:,iel)-1)+2;
    ind_P            = (np*(iel-1)+1):(np*(iel-1)+np);
    
    A_all(ind_el, ind_el) = A_all(ind_el, ind_el) + A_elem;
    Q_all(ind_P, ind_el)  = Q_all(ind_P, ind_el) + Q_elem';
    invM_all(ind_P,ind_P) = invM_all(ind_P,ind_P) + invM_elem;
    Rhs_all(ind_el)       = Rhs_all(ind_el) + Rhs_elem(:);
end

A = A_all;
Q = Q_all;
invM = invM_all;
Rhs = Rhs_all;

Vel            = zeros(2*nnod,1);
Vel(bc_ind)    = bc_val;

Free        = 1:sdof;
Free(bc_ind)= [];
Rhs         = Ra*Rhs -  A*Vel;
A           = A(Free,Free);

L = chol(A, 'lower');
Pressure    = zeros(nel*np, 1);
while (div_max>div_max_uz  && uz_iter<uz_iter_max)
    uz_iter         = uz_iter + 1;
    Vel(Free)       = L'\(L\(Rhs(Free)));
    %     Vel(Free)         = A\Rhs(Free);
    Div             = invM*(Q*Vel);                                        %COMPUTE QUASI-DIVERGENCE
    Rhs             = Rhs - PF*(Q'*Div);                                   %UPDATE RHS
    Pressure        = Pressure + PF*Div;                                   %UPDATE TOTAL PRESSURE (negative sign convention)
    div_max         = max(abs(Div(:)));                                    %CHECK INCOMPRESSIBILITY
end

test = 0;















