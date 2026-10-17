function [num_components,N] = extract_connected_components_bipartite(N,Spt0)
    [m, n] = size(N);    
    component_ids=findcomp(N);
    unique_components = unique(component_ids);
    num_components = length(unique_components);
    if num_components==1
        return
    else
    com = find(component_ids>1);
    rows=com(com<m+1);
    cols=com(com>m)-m;
    row_mask = ismember(1:size(Spt0,1), rows);
    col_mask = ismember(1:size(Spt0,2), cols);

    N(row_mask, :) = N(row_mask, :) | Spt0(row_mask, :);
    N(:, col_mask) = N(:, col_mask) | Spt0(:, col_mask);
    component_ids=findcomp(N);
    unique_components = unique(component_ids);
    num_components = length(unique_components);
    end
end 
function component_ids=findcomp(N)
    [m, n] = size(N);
    total_nodes = m + n;
    A = sparse(total_nodes, total_nodes);
    A(1:m, m+1:end) = N;         
    A(m+1:end, 1:m) = N';        
    G = graph(A, 'upper', 'OmitSelfLoops');
    component_ids = conncomp(G, 'OutputForm', 'vector');
end