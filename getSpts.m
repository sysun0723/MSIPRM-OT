function Spts = getSpts(Spt, lambda1, lambda2)
    N=size(Spt,1);
    n=single(sqrt(N));
    lambda1 = single(lambda1(:));
    lambda2 = single(lambda2(:));
    N = single(n*n);
    if n < 129 
        block_size = single(1024);
    elseif 255<n  &&  n<513
        block_size = single(4096);
    elseif n==1024
        block_size = single(8192);
    elseif n==2048
        block_size = single(16384);
    else
        block_size = single(4096);
    end 
    
    num_blocks = ceil(N/block_size);
    Spt_i = cell(num_blocks, num_blocks);
    Spt_j = cell(num_blocks, num_blocks);
    [rows, cols] = size(Spt);
    row_blocks = ceil(rows / block_size);
    col_blocks = ceil(cols / block_size);
    [i,j] = find(Spt);
    block_i = floor((i-1)/block_size) + 1;
    block_j = floor((j-1)/block_size) + 1;
    block_map = accumarray([block_i, block_j], 1, [row_blocks, col_blocks], @any, false);
    [blocks(:,1), blocks(:,2)] = find(block_map);
    
    for idx = 1:size(blocks, 1)
        k = blocks(idx, 1);
        m = blocks(idx, 2);
        
        i_range = single((k-1)*block_size+1 : min(k*block_size, N))';    
        j_range = single((m-1)*block_size+1 : min(m*block_size, N)) ;
           
        i1 = mod(i_range-1, single(n));
        i2 = floor((i_range-1)/single(n));
        j1 = mod(j_range-1, single(n));
        j2 = floor((j_range-1)/single(n));
        vv = single((2.*i1.*j1 + 2.*i2.*j2));

            
        ii = (i1.^2 + i2.^2) - lambda2(i_range) ;
        jj = (j1.^2 + j2.^2) - lambda1(j_range)';
        uu = single(ii + jj);
        
        %mask=uu<vv;
        [i,j] = find(uu<vv);
        Spt_i{k,m} = i_range(i) ;
        Spt_j{k,m} = j_range(j)';
    end
    Spts = sparse(double(vertcat(Spt_i{:})), double(vertcat(Spt_j{:})), true, double(N), double(N));
end
