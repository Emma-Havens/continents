function T = thermal2d_transient(ELEM2NODE, PHASES, GCOORD, D, Bc_ind, nip, T, dt, elemType, reorder, method, lumping)
% THERMAL2D Two dimensional finite element thermal problem solver of MILAMIN

%   Part of MILAMIN: MATLAB-based FEM solver for large problems, Version
%   1.0.1
%   Copyright (C) 2007, M. Dabrowski, M. Krotkiewski, D.W. Schmid
%   University of Oslo, Physics of Geological Processes
%   http://milamin.org
%   See License file for terms of use.

%==========================================================================
% MODEL INFO
%==========================================================================
nnod         = size(GCOORD,2);
nnodel       = size(ELEM2NODE,1);
nel          = size(ELEM2NODE,2);

%==========================================================================
% CONSTANTS
%==========================================================================
ndim         =   2;
nelblo       = 760;

%==========================================================================
% BLOCKING PARAMETERS (nelblo must be < nel)
%==========================================================================
nelblo       = min(nel, nelblo);
nblo         = ceil(nel/nelblo);

%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
switch elemType
    case 'quad'
        %         [IP_X, IP_w] = ip_quad(nip);
        [IP_X, IP_w] = ip_line(nip);
        [xx yy] = meshgrid(IP_X,IP_X);
        IP_X = [xx(:) yy(:)];
        IP_w = IP_w'*IP_w;
        IP_w = IP_w(:);
        nip  = length(IP_w);
        [N     dNdu] = shp_deriv_quad(IP_X,nnodel);
    case 'triangle'
        [IP_X, IP_w] = ip_triangle(nip);
        [N     dNdu] = shp_deriv_triangle(IP_X,nnodel);
end

%==========================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==========================================================================
K_all        = zeros(nnodel*(nnodel+1)/2,nel); 
M_all        = zeros(nnodel*(nnodel+1)/2,nel);

%==========================================================================
% INDICES EXTRACTING LOWER PART
%==========================================================================
indx_l       = tril(ones(nnodel)); indx_l = indx_l(:); indx_l = indx_l==1;

%==========================================================================
% CALCULATING THE MASS MATRIX
%==========================================================================


switch method
    %======================================================================
    % STANDARD VERSION
    %======================================================================        
    case 'std'        
        %==================================================================
        % DECLARE VARIABLES (ALLOCATE MEMORY)
        %==================================================================
        K_elem      = zeros(nnodel,nnodel);   
        M_elem      = zeros(nnodel,nnodel);
        
        %==================================================================
        % i) ELEMENT LOOP - MATRIX COMPUTATION
        %==================================================================
        for iel = 1:nel
            %==============================================================
            % ii) FETCH DATA OF ELEMENT
            %==============================================================
            ECOORD_X = GCOORD(:,ELEM2NODE(:,iel));
            ED       = D(PHASES(iel));

            %==============================================================
            % iii) INTEGRATION LOOP
            %==============================================================
            K_elem(:) = 0;
            M_elem(:) = 0;
            for ip=1:nip
                %==========================================================
                % iv) LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
                %==========================================================
                Ni          = N{ip};
                dNdui       = dNdu{ip};

                %==========================================================
                % v) CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
                %==========================================================
                J           = ECOORD_X*dNdui;
                detJ        = det(J);
                invJ        = inv(J);

                %==========================================================
                % vi) DERIVATIVES wrt GLOBAL COORDINATES
                %==========================================================
                dNdX        = dNdui*invJ;

                %==========================================================
                % vii) NUMERICAL INTEGRATION OF ELEMENT MATRICES
                %==========================================================
                K_elem      = K_elem + IP_w(ip)*detJ*ED*(dNdX*dNdX');
                M_elem      = M_elem + IP_w(ip)*detJ*(Ni*Ni');
            end

            %==============================================================
            % ix) WRITE DATA INTO GLOBAL STORAGE
            %==============================================================
            K_all(:,iel)    = K_elem(indx_l);
            M_all(:,iel)    = M_elem(indx_l);
        end
                
    %======================================================================
    % OPTIMIZED VERSION
    %======================================================================    
    case 'opt'        
        %==================================================================
        % DECLARE VARIABLES (ALLOCATE MEMORY)
        %==================================================================        
        K_block     = zeros(nelblo,nnodel*(nnodel+1)/2);
        M_block     = zeros(nelblo,nnodel*(nnodel+1)/2);
        invJx       = zeros(nelblo, ndim);
        invJy       = zeros(nelblo, ndim);
        il          = 1;
        iu          = nelblo;

        %==================================================================
        % i) BLOCK LOOP - MATRIX COMPUTATION
        %==================================================================
        for ib = 1:nblo
            %==============================================================
            % ii) FETCH DATA OF ELEMENTS IN BLOCK
            %==============================================================
            ECOORD_x = reshape( GCOORD(1,ELEM2NODE(:,il:iu)), nnodel, nelblo);
            ECOORD_y = reshape( GCOORD(2,ELEM2NODE(:,il:iu)), nnodel, nelblo);
            ED       = reshape(D(PHASES(il:iu)),nelblo,1);

            %==============================================================
            % iii) INTEGRATION LOOP
            %==============================================================
            K_block(:)  = 0;
            M_block(:)  = 0;
            for ip=1:nip
                %==========================================================
                % iv) LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
                %==========================================================
                Ni          = N{ip};
                dNdui       = dNdu{ip};

                %==========================================================
                % v) CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
                %==========================================================
                Jx          = ECOORD_x'*dNdui;
                Jy          = ECOORD_y'*dNdui;
                detJ        = Jx(:,1).*Jy(:,2) - Jx(:,2).*Jy(:,1);

                invdetJ     = 1.0./detJ;
                invJx(:,1)  = +Jy(:,2).*invdetJ;
                invJx(:,2)  = -Jy(:,1).*invdetJ;
                invJy(:,1)  = -Jx(:,2).*invdetJ;
                invJy(:,2)  = +Jx(:,1).*invdetJ;

                %==========================================================
                % vi) DERIVATIVES wrt GLOBAL COORDINATES
                %==========================================================
                dNdx        = invJx*dNdui';
                dNdy        = invJy*dNdui';

                %==========================================================
                % vii) NUMERICAL INTEGRATION OF ELEMENT MATRICES 
                %==========================================================
                weight      = IP_w(ip)*detJ;

                indx = 1;
                for i = 1:nnodel
                    for j = i:nnodel
                        K_block(:,indx)  =   K_block(:,indx) + ...
                            (dNdx(:,i).*dNdx(:,j)+ dNdy(:,i).*dNdy(:,j)).*weight.*ED;
                        M_block(:,indx)  =   M_block(:,indx) + ...
                            (Ni(i).*Ni(j)).*weight;
                        indx = indx + 1;
                    end
                end
            end
            %==============================================================
            % ix) WRITE DATA INTO GLOBAL STORAGE
            %==============================================================
            K_all(:,il:iu)	= K_block';
            M_all(:,il:iu)	= M_block';
            
            %==============================================================
            % READJUST START, END AND SIZE OF BLOCK. REALLOCATE MEMORY
            %==============================================================
            il  = il+nelblo;
            if(ib==nblo-1)
                nelblo 	= nel-iu;
                K_block	= zeros(nelblo, nnodel*(nnodel+1)/2);
                M_block	= zeros(nelblo, nnodel*(nnodel+1)/2);
                invJx   = zeros(nelblo, ndim);
                invJy   = zeros(nelblo, ndim);
            end
            iu  = iu+nelblo;
        end
end

%==========================================================================
% ix) CREATE TRIPLET FORMAT INDICES
%==========================================================================
indx_j = repmat(1:nnodel,nnodel,1); indx_i = indx_j';
indx_i = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
indx_j = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);

K_i = ELEM2NODE(indx_i,:); K_i = K_i(:); 
K_j = ELEM2NODE(indx_j,:); K_j = K_j(:);

indx       = K_i < K_j;
tmp        = K_j(indx);
K_j(indx)  = K_i(indx);
K_i(indx)  = tmp;

%==========================================================================
% x) CONVERT TRIPLET DATA TO SPARSE MATRIX
%==========================================================================
K_all  = K_all(:);
M_all  = M_all(:);
if exist(['sparse2.' mexext], 'file') == 3
    K      = sparse2(K_i, K_j, K_all);
    M      = sparse2(K_i, K_j, M_all);
else
    K      = sparse(double(K_i), double(K_j), K_all);
    M      = sparse(double(K_i), double(K_j), M_all);
end
clear K_i K_j K_all;

%==========================================================================
% LUMPING
%==========================================================================
if lumping
    M  = M + M' - spdiags(spdiags(M,0),0,length(M),length(M));
    M  = spdiags(sum(M,2), 0,length(M),length(M));     
end

% Fully implicit
K      = dt*K + M;

clear K_i K_j K_all M_all;

%==========================================================================
% BOUNDARY CONDITIONS
%==========================================================================
if exist(['cs_transpose.' mexext], 'file') == 3
    TMP         = K(:,Bc_ind) + cs_transpose(K(Bc_ind,:));
else
    TMP         = K(:,Bc_ind) + K(Bc_ind,:)';
end
Rhs         = M*T + M'*T - spdiags(M,0).*T;
Rhs         = Rhs - TMP*T(Bc_ind);  

Free         = 1:nnod;
Free(Bc_ind) = [];
K            = K(Free,Free);

%==========================================================================
% REORDERING
%==========================================================================
switch reorder
    case 'metis'
        perm = metis(K);
    case 'amd'
        perm = amd(K);
    otherwise
        error('Unknown reordering')
end

%==========================================================================
% FACTORIZATION - ideally L = lchol(K, perm)
%==========================================================================
if exist(['cs_transpose.' mexext], 'file') == 3 && exist(['cs_symperm.' mexext], 'file') == 3
    K = cs_transpose(K);
    K = cs_symperm(K,perm);
    K = cs_transpose(K);
    L = lchol(K);
else
    K = K(perm,perm);
    K = tril(K)+triu(K,1)';
    L = chol(K, 'lower');
end


%==========================================================================
% SOLVE
%==========================================================================
if exist(['cs_ltsolve.' mexext], 'file') == 3 && exist(['cs_lsolve.' mexext], 'file') == 3
    T(Free(perm)) = cs_ltsolve(L,cs_lsolve(L,Rhs(Free(perm))));
else
    T(Free(perm)) = L'\(L\(Rhs(Free(perm))));
end

