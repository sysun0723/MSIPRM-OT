function [Spt,Spt_i, Spt_j, R, p] = enlarge_Spt(Z, m, n, rows, cols)
    Spt0=sparse(Z>0);
    N=nnz(Spt0);
    [S_i,S_j]=find(Spt0);
    k0=2;
    z = Z(Spt0);
    [~, idx] = sort(z, 'descend');   
    k = min(round(k0 * (m + n)), N);
    zz = zeros(size(z));     
    zz(idx(1:k)) = z(idx(1:k)); 
    Z=sparse(S_i,S_j,zz,size(Spt0,1),size(Spt0,2)); 
    Spt = Z >0;
    [num_components,Spt] = extract_connected_components_bipartite(Spt,Spt0); 
    Spt_expanded = expandSupport(Spt, rows/2, cols/2);
    [Spt_i, Spt_j] = find(Spt_expanded);
    [R, p] = Permutation_matrix(Spt_expanded, Spt_i, Spt_j);
    Spt = Spt_expanded;
end

function Spt = expandSupport(Spt, m, n) 
    newSize = 4*m*n;
    [Spt_i,Spt_j] = find(Spt); 
    i=4*Spt_i-2*(mod(Spt_i - 1, m) + 1)-1;
    j=4*Spt_j-2*(mod(Spt_j - 1, m) + 1)-1;
    u=[i+0;i+1;i+2*m;i+2*m+1];
    v=repmat(j,4,1);
    uu=repmat(u,4,1);
    vv=[v+0;v+1;v+2*m;v+2*m+1];
    ij=[uu,vv];
    values = true(size(ij,1),1);  
    Spt = sparse(ij(:,1), ij(:,2), values, newSize, newSize); 
end

function [R,p]=Permutation_matrix(Spt,Spt_i,Spt_j)
        V1 = sum(Spt,1);   
        V2 = sum(Spt,2);
        V2 = V2(1:end-1); 
        v=(V1(Spt_j))';
        dd=1./sqrt(v);
        DD= sparse(Spt_i, Spt_j, dd, size(Spt,1), size(Spt,2));
        DD0=DD(1:end-1, :);
        main_func = diag(V2) - DD0 * DD0';
        [R,~,p] = chol(main_func,'vector');
end

