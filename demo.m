clc; 
clear;
%matrix=loadData(type,size,num);
% type=1: CauchyDensity;   type=2: ClassicImages;   type=3: GRFmoderate;        type=4: GRFrough;    type=5: GRFsmooth
% type=6: LogGRF;          type=7: LogitGRF;        type=8: MicroscopyImages;   type=9: Shapes;      type=10: WhiteNoise
%size=32,64,128,256,512;  num=1001,...,1010; type=1,...,10
% loadData(1,32,1003) reads DOTmark_1.0/Data/CauchyDensity/data32_1003.csv beside this script.

type= 2 ;
num1=1004;
num2=1005;
% num_start = 1001;
% num_end = 1010;
% pairs = nchoosek(num_start:num_end, 2);
% pairs = [1002 1005; 1001 1007; 1003 1008; 1004 1010; 1010 1009;  
%          1001 1003; 1005 1010; 1002 1008; 1004 1007; 1006 1009];
% pairs = [1002 1010; 1001 1007; 1007 1008; 1004 1005; 1006 1009];
% 
% all_results = cell(10, 1);  
% 
% for type = [5]
%     fprintf('Processing type=%d ...\n', type);
%     result = zeros(size(pairs,1), 10);   
% for tt = 1:size(pairs,1)
%     num1 = pairs(tt,1); 
%     num2 = pairs(tt,2);

scale = 512;
matrix_1= loadData(type,scale, num1);
matrix_2= loadData(type,scale, num2);
maxit = 200;
printlevel=1;
% scale=1024;
multiscale_max=scale/16;
k=multiscale_max;
elapsed_time=tic;
    matrix_1_1=matrix_scale(matrix_1,k);
    matrix_2_1=matrix_scale(matrix_2,k);
    [m,n,Su,De]=basic_data(matrix_1_1,matrix_2_1);
    b = [Su; De];
    [x,s,lambda,z,opt,iter] = iprm4denseot(m,n,Su,De,printlevel);
    z(z<1e-13)=0;
    Z=reshape(z,[m,n]);
    k=k/2;
    Spt=sparse(Z>0);
while k>=1 
    matrix_1_1=matrix_scale(matrix_1,k);
    matrix_2_1=matrix_scale(matrix_2,k);
    [rows,cols]=size(matrix_1_1);
    % if rows<65  && rows<scale
    % tol1=1e-6; 
    % elseif 64<rows  && rows<257
    % tol1=1e-6;   
    % elseif rows==256
    % tol1=1e-6; 
    % elseif rows==512
    %     tol1=1e-6;
    %     tol=1e-6;
    % else
    % tol1=1e-6;     
    % end
     
    [Spt, Spt_i, Spt_j,R, p] = enlarge_Spt(Z, m, n, rows, cols);  
    norm_c=normc(rows);
    m=rows*cols;
    n=rows*cols;
    Su=matrix_1_1(:);
    De=matrix_2_1(:);
    Su = Su/sum(Su);
    De = De/sum(De);
    % nzSu=Su>0;
    % nzDe=De>0;
    % Su=Su(nzSu);
    % De=De(nzDe);
    b = [Su; De];
    b=b/norm(b);

    c_Spt = construct_c_Spt(rows,cols,Spt_i,Spt_j);
  
    [x,s,lambda]= start_Spt(m,n,b,c_Spt,R,p,Spt,Spt_i,Spt_j);
 
    update_Spt_num=1;
    usestart=1;
    mu=1;
    rho=1;
    while update_Spt_num<10
    if rows<scale/2+1
       tol1=1e-6;
       tol=1e-6; 
    elseif rows==scale
       tol1=1e-6;
       tol=1e-6; 
    end

    [Spt_i,Spt_j]=find(Spt);
    c_Spt=construct_c_Spt(rows,cols,Spt_i,Spt_j); 
    if  usestart==0
    x=X(Spt);
    s= c_Spt-(lambda1(Spt_j)+lambda2(Spt_i));
    [R,p]=sparseChol(Spt,Spt_i,Spt_j);
    end
    usestart=0;
    
   % [connected_components] = extract_connected_components_bipartite(Spt)
    nnz(Spt)

   solvestart=tic;
    [x,s,lambda,z,residue,iter0,mu,rho,opt] = iprm4sparseot(m,n,c_Spt,b,x,s,lambda,R,p,Spt,mu,rho,tol1,maxit,printlevel);
   solve_time = toc(solvestart);  
   disp(['solvetime=', num2str(solve_time), 's']); 

    iter=iter+iter0;
    X = sparse(Spt_i, Spt_j, x, size(Spt,1), size(Spt,2));
    S = sparse(Spt_i, Spt_j, s, size(Spt,1), size(Spt,2));
    Z = sparse(Spt_i, Spt_j, z, size(Spt,1), size(Spt,2));
       Spt0= Z>0;
       lambda=[lambda;0];          %#ok 
       lambda1=lambda(1:m);
       lambda2=lambda(m+1:end);

       getstart=tic;
       Spts=getSpts(Spt,lambda1,lambda2);
       getspts_time=toc(getstart);
       disp(['getsptstime=', num2str(getspts_time), 's']);


       Spts=Spts-Spt>0;
       [Spts_i,Spts_j]=find(Spts);

       ll=lambda1(Spts_j)+lambda2(Spts_i);
       c_Spts=construct_c_Spt(rows,cols,Spts_i,Spts_j);
       s0= c_Spts-ll;
       residue(3)=sqrt(s0'*s0+(residue(3)*(1+norm(c_Spt)))^2)/(1+norm_c);
       %disp(residue(3))

    if residue(3)<tol   
        disp('problem solved');
        break;
    else
        tol1=max(1e-1*tol1,1e-8);
        %C_Spts=sparse(Spts_i, Spts_j, c_Spts, size(Spts,1), size(Spts,2));
        %Spts=id_Spts(Spts,Spts_i,Spts_j, c_Spts,ll,nnz0,ids_num);
        ids_num=0.2;
     if nnz(Spts)>2*ids_num*(m+n) 
      s=c_Spts./ll;
      l=ceil(2*ids_num*(m+n));
      [~, idx] = mink(s, l);       
      binary_vec = zeros(size(s)); 
      binary_vec(idx) = 1;
      Spts=sparse(Spts_i, Spts_j, binary_vec, size(Spts,1), size(Spts,2));
     end
     Spt=Spts+Spt0>0;
    end
    
     pinf =norm(min(x,0))/(1+norm(x));
     gap = residue(2);
     dinf = residue(3);
    
     disp(opt(2));
     disp(opt(3));
     disp(pinf);
     disp(gap );
     disp(dinf) 


    update_Spt_num=update_Spt_num+1;
    end
    k=k/2;
end
total_time=toc(elapsed_time);
disp(['total_time=', num2str(total_time), 's']);

 %pinf = residue(1);
 gap = residue(2);
 dinf = residue(3);

 disp(opt(2));
 disp(opt(3));
 disp(pinf);
 disp(gap );
 disp(dinf) 

% fprintf('Finished %d/%d: (%d,%d), Solvetime %.2f秒\n',...
%             tt, size(pairs,1), num1, num2, total_time);
% result(tt,:) = [num1, num2, total_time, opt(1), iter,opt(2),opt(3), pinf, gap, dinf]; 
% end
% 
% result = array2table(result, ...
%     'VariableNames', {'num1', 'num2', 'time', 'opt','iter','<c,x>', '<b,lambda>', 'pinf', 'gap', 'dinf'});
% disp(result);
% 
%     save_path = fullfile(fileparts(mfilename('fullpath')), 'data', 'MSIPRM', '1024');
%     filename = fullfile(save_path, sprintf('matlab-%d-%d.csv',  type, scale));
%     writetable(result, filename);
%     fprintf('type=%d The result has been saved to‌ %s\n', type, filename);
% end 


function importedData = loadData(type, size, num)
    basePath = fullfile(fileparts(mfilename('fullpath')), 'DOTmark_1.0', 'Data');
    folderMap = containers.Map( [1, 2, 3,4,5,6,7,8,9,10], {'CauchyDensity', 'ClassicImages', 'GRFmoderate',...
        'GRFrough','GRFsmooth','LogGRF','LogitGRF','MicroscopyImages','Shapes','WhiteNoise'});
    folderName = folderMap(type);
    filename = fullfile(basePath, folderName, sprintf('data%d_%d.csv',  size, num));
    importedData = readmatrix(filename);
end 

function C = construct_C(m,n)
    total = m * n;                      
    indices = (1:total)';               
    L = mod(indices - 1, n) + 1;       
    K = floor((indices - 1) / m) + 1;  
    row_diff = L - L';                  
    col_diff = K - K';                 
    C = row_diff.^2 + col_diff.^2;
end 

function c = construct_c_Spt(m,n,Spt_i,Spt_j)       
  c=(mod(Spt_i-1,n)-mod(Spt_j-1,n)).^2+(floor((Spt_i-1)/m)-floor((Spt_j-1)/m)).^2;
end 

function [m,n,Su,De]=basic_data(matrix_1,matrix_2)
[rows,cols]=size(matrix_1);
m=rows*cols;
n=rows*cols;
Su=matrix_1(:);
De=matrix_2(:);
Su = Su/sum(Su);
De = De/sum(De); 
end

% function X= northwest(m,n,Su,De) 
%     row_indices = [];
%     col_indices = [];
%     values = [];
%     i = 1; j = 1;
%     while i <= m && j <= n 
%         alloc = min(Su(i), De(j));
%         if alloc > 0   
%             row_indices(end+1) = i;
%             col_indices(end+1) = j;
%             values(end+1) = alloc;
%         end 
%         Su(i) = Su(i) - alloc;
%         De(j) = De(j) - alloc;
%         if Su(i) == 0
%             i = i + 1;
%         else 
%             j = j + 1;
%         end
%     end 
%     X = sparse(row_indices, col_indices, values, m, n);
%     X=X';
% end


function [x,s,lambda]= start_Spt(m,n,b,c,R,p,Spt,Spt_i,Spt_j) 
        V1 = sum(Spt,1)';   
        St=Spt(1:end-1, :);
        cc=matvec_spt(c,Spt, Spt_i,Spt_j);
        c1 = cc(1:m);
        c2 = cc(m+1:m+n-1);
        rhsSchur = c2 - St*(c1./V1);
        la2=zeros(size(rhsSchur));
        la2(p) = R\(R'\(rhsSchur(p)));
        la1 = (c1 - St'*la2) ./ V1;
        lambda= [la1; la2];
        s=c-matTvec_spt(m,n,lambda, Spt_i,Spt_j);
        b1=b(1:m);
        b2=b(m+1:end-1);
        rhsSchur = b2 - St*(b1./V1);
        la2=zeros(size(rhsSchur));
        la2(p) = R\(R'\(rhsSchur(p)));
        la1 = (b1 - St'*la2) ./ V1;
        la= [la1; la2];
        x=matTvec_spt(m,n,la,Spt_i,Spt_j);
end

function [R,p]=sparseChol(Spt,Spt_i,Spt_j)
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

  function result = matvec_spt(x,Spt, Spt_i,Spt_j)
  X = sparse(Spt_i, Spt_j, x, size(Spt,1), size(Spt,2));
  x1=sum(X,1)';
  x2=sum(X,2);
  x2=x2(1:end-1);
  result=[x1;x2];
  end

function result = matTvec_spt(input_m, input_n, lambda, Spt_i,Spt_j)       
  lambda=[lambda;0];
  lambda1=lambda(1:input_m);
  lambda2=lambda(input_n+1:end);
  result=lambda1(Spt_j)+lambda2(Spt_i);
end 

function result=matrix_scale(A,k)
[m, n] = size(A); 
if k==1
    result=A;
else
new_m = m / k; new_n = n / k;
B = reshape(A, [k, new_m, k, new_n]);
result= squeeze(mean(mean(B, 1), 3));
end
end

function norm_c = normc(m)
    norm_c=m^2*sqrt((m^2-1)*(17*m^2-23)/90);
end

function A = construct_A_Spt(n, m,Spt)
     jj=find(Spt>0);
     ii1=floor((jj-1)./m)+1;
     ii2=m+mod(jj-1,m)+1;
     jj= (1:nnz(Spt))';
     A = sparse([ii1; ii2], [jj; jj],true,m+n,nnz(Spt));
     A = A(1:end-1, :);
end
