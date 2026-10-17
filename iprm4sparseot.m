function [x,s,lambda,z,residue,iter,mu,rho,opt] = iprm4ot2(input_m,input_n,c_Spt,b,x,s,lambda,R,p,Spt,mu,rho,tol,maxit,printlevel)
b0=b(1:end-1);
[Spt_i,Spt_j]=find(Spt);
b1=b(1:input_m);
b2=b(input_m+1:end);
if (size(b0,2) > 1) 
    b0 = (b0)';
end
if (issparse(b0))  
    b0 = full(b0);  
end

pl = printlevel;
opt=0;
gamma0=1e-3;
epsi=tol;    %stopping parameter
delta=0.5;    %stepsize parameter
sigma=0.01;   %line search parameter
eta = 5;
beta=0.5;     %update rho parameter      0.05  0.1 0.3  0.5

if input_m<32^2+1
    usepcg =1;   %use pcg
    usechol=0;  %use chol
else
    usepcg = 0;   %use pcg
    usechol=1;  %use chol
end
mu=0.1;
rho=10;
[x,s,lambda]= Warmstart_Spt(input_m,input_n,b,c_Spt,x,lambda,s,mu,rho,R,p,Spt,Spt_i,Spt_j);

[z,t] = compz(rho,mu,x,s);
%y=mu./z;
res=[matTvec_Spt(input_m,input_n,lambda, Spt_i,Spt_j)+s-c_Spt;z-x;matvec_Spt0(x,Spt,Spt_i,Spt_j)-b];
% res=[matTvec(input_m,input_n,lambda)+s-c;z(pos_vars)-x(pos_vars);matvec(input_m,input_n,x)-b];
phi=0.5*(res'*res);
norm_c=norm(c_Spt);

while mu > eta*max(epsi/eta,phi)
    mu = mu /eta;
    [z,t] = compz(rho,mu,x,s);
    %res=[matTvec(input_m,input_n,lambda)+s-c;z-x;matvec(input_m,input_n,x)-b]; 
    %res=[z-x;matvec(input_m,input_n,x)-b]; 
    res=z-x; 
    phi=0.5*(res'*res);
end
iter=0;
header(pl);
while iter<maxit 
    if  7<iter 
      rh=zeros(size(b));
      rd=zeros(size(s));
    else
    rh = matvec_Spt0(x,Spt,Spt_i,Spt_j)-b;
    rd =c_Spt-matTvec_Spt(input_m,input_n,lambda, Spt_i,Spt_j)-s;
    end
    norm_rh=norm(rh);
    norm_rd=norm(rd);
%     norm_res=[norm(rh)/(1+norm_b);norm(rd)/(1+norm_c);abs(c'*x-b0'*lambda)/(1+abs(c'*x)+abs(b0'*lambda))];
%     norm_res=[norm(rh)/(1+norm_b);norm(rd)/(1+norm_c);mu];
if input_m<64^2+1
    norm_res=[abs(c_Spt'*x-b0'*lambda)/(1+abs(c_Spt'*x)+abs(b0'*lambda));norm_rh/2;norm_rd/(1+norm_c);100*mu;100*phi];
elseif input_m>64^2  &&  input_m<128^2+1
    norm_res=[abs(c_Spt'*x-b0'*lambda)/(1+abs(c_Spt'*x)+abs(b0'*lambda));norm_rh/2;norm_rd/(1+norm_c);10*mu;10*phi];
elseif input_m>=256^2
    norm_res=[abs(c_Spt'*x-b0'*lambda)/(1+abs(c_Spt'*x)+abs(b0'*lambda));norm_rh/2;norm_rd/(1+norm_c);mu;phi];
end
    % if (nnz(isnan(rh)) > 0 || nnz(isinf(rh)) > 0 || nnz(isnan(rd)) > 0 || nnz(isinf(rd)) > 0)
    %     opt = -1;
    %     return;
    % end
    if isempty(find(norm_res >= epsi,1)) 
        opt=1;
        disp('problem solved');
%        [norm(rh)/(1+norm_b);norm(rd)/(1+norm_c);abs(c'*x-b0'*lambda)/(1+abs(c'*x))];
        break;
    end    
    
    %The search direction
%     norm_rb=norm(rh);
%     norm_rc=norm(rd);
%  if mu>1e-6
%   du =-mu +0.5*gamma*phi;
% elseif mu<1e-6 || phi<1e-7
% du =-mu +0.5*gamma*phi;
% % elseif 1e3*norm_rb<mu || 1e3*norm_rc <mu
% % du=-mu +0.75*gamma*phi;
% % elseif 1e2*norm_rb<mu || 1e2*norm_rc<mu
% % du =-mu +0.75*gamma*phi; 
% else
%  du=-mu+0.995*gamma*phi;
% end    
if  mu<1e-5  && phi>1e+0
    gamma0=1e-3;
elseif    mu>3e-4  &&  iter<6
    gamma0=5e-3;
elseif mu>3e-4 && iter>=5
    gamma0=1e-3;
elseif mu>3e-5 && iter>=5  &&  mu<3e-4
    gamma0=5e-4;
elseif  mu<3e-5  &&  mu>3e-6
    gamma0=2e-4;
elseif  mu<3e-6  &&  phi>3e-6
    gamma0=1e-3;
elseif  mu<3e-6  &&  phi<3e-6  && phi>1e-9
    gamma0=1e-2;
% elseif  mu<3e-6  &&  mu>1e-8   &&  phi<1e-9
%     gamma0=3e-4;
%     du=-mu+gamma*phi;
% elseif  mu<1e-8  
%     gamma0=2e-5;
%     du=-mu+gamma*phi;
end
    gamma = min(gamma0,mu/phi);
    du=-mu+gamma*phi;
    rf = du+(t.*(z-x));
    
    if  5<iter && iter<80
    rhs = mu*x+z.*rf;
   % rhs = mu*b0-matvec_Spt(mu*x+z.*rf,Spt,Spt_i,Spt_j);
    else
    rhs = mu*x+z.*(-z.*rd+rf);
    end
    Rhs = sparse(Spt_i, Spt_j, rhs, size(Spt,1), size(Spt,2));
    rhs1= mu*b1-sum(Rhs,1)';
    rhs2= mu*b2-sum(Rhs,2);
    %zz = z.*z;
    ZZ = sparse(Spt_i, Spt_j, z.*z, size(Spt,1), size(Spt,2));
    V1 = (sum(ZZ,1));
    V2 = (sum(ZZ,2));

    % rhs1 = rhs(1:input_m);
    % rhs2 = rhs(input_m+1:input_m+input_n-1);
    %rhsSchur = rhs2 - ZZ*(rhs1./V1');

if usechol 
        %ZZ(end, :) = [];
        DD = ZZ(p,:);
        D=DD./sqrt(V1);
        rhsSchur = rhs2(p) - DD*(rhs1./V1');

        main_func= diag(V2(p)) - D*D';
        % [R,flag] = chol(main_func);
        % if flag~=0
        %    opt=0;
        %    break
        % end
        % dla2=zeros(size(rhsSchur));
        % dla2(p)=(R\(R'\(rhsSchur(p)))); 
        dla2=zeros(size(rhsSchur));
        dla2(p) = main_func\(rhsSchur); 
        dla1 = (rhs1 - ZZ'*[dla2;0]) ./ V1';
end 

if usepcg 
        ZZ(end, :) = [];
        rhs2= rhs2(1:end-1);
        rhsSchur = rhs2 - ZZ*(rhs1./V1');
        V2 = V2(1:end-1);
        main_func2 = @(x)MAIN_FUNC2(ZZ, V1, V2, x); 
        pre_fun2 = @(x)PRE_FUN2(V2, x);
        maxitpcg = 500;
        pcgtol = min(1e-8, 1e-4*norm(rhsSchur));
        [dla2, flag,~, iter_pcg] = pcg(main_func2, rhsSchur, pcgtol, maxitpcg, pre_fun2);
         dla1 = (rhs1 - ZZ'*dla2) ./ V1';
        if iter_pcg >500  || flag~=0
            usepcg=0;
            usechol=1;
        end

        % DD = ZZ(p,:);
        % D=DD./sqrt(V1);
        % rhsSchur = rhs2(p) - DD*(rhs1./V1');
        % 
        % main_func= diag(V2(p)) - D*D';
        % opts.type = 'ict';           
        % opts.droptol = 1e-3;
        % L = ichol(main_func, opts);
        % pcgtol = min(1e-8, 1e-4*norm(rhsSchur));
        % dla2=zeros(size(rhsSchur));
        % main_func2=@(x)MAIN_FUNC(D, V2,p, x);
        %  [dla2(p), flag,~, iter_pcg] = pcg(main_func2, rhsSchur, pcgtol, maxitpcg, L,L');
        % 

        %[dla2(p), flag,~, iter_pcg] = pcg(main_func, rhsSchur, pcgtol,
        %maxitpcg, L,L');
        %dla1 = (rhs1 - ZZ'*[dla2;0]) ./ V1';

        if iter_pcg >500  || flag~=0
            usepcg=0;
            usechol=1;
        end


end


% if useldl 
%         sqV1=sqrt(V1);
%         v=(sqV1(Spt_j))';
%         d=zz./v;
%         D= sparse(Spt_i, Spt_j, d, input_m, input_n);
%         D=D(p,:);
%         main_func= diag(V2(p)) - D*D';
%         [R,flag] = chol(main_func);
%          [L, D, P] = ldl(main_func,'vector');
%         if flag~=0
%            opt=0;
%            break
%         end
%         dla2=zeros(size(rhsSchur));
%         dla2(p)=(R\(R'\(rhsSchur(p))));  
% end 

    dla = [dla1; dla2];
    
    % if iter>7 
    % ds = -matTvec_Spt(input_m,input_n,dla, Spt_i,Spt_j);
    % else
    ds = rd-matTvec_Spt(input_m,input_n,dla, Spt_i,Spt_j);
   % end
    dx = (z.*(rf-z.*ds))./mu;

    %Doing line search
if mu>3e-4  &&  iter<6
    alpha=1;
elseif mu>3e-4 && iter>=5 
    alpha=0.25;
elseif mu>3e-5 && iter>=5  &&  mu<3e-4
    alpha=0.25;
elseif  mu<3e-6  &&  phi>3e-7
    alpha=0.5;
elseif  mu<3e-6  &&  phi<3e-7 
    alpha=1;
else
    alpha=1;
end

% if mu>3e-4  &&  iter<10
%     alpha=1;
% elseif mu>3e-4 && iter>=10 
%     alpha0=0.25;
%     if alpha<=alpha0/4
%         alpha=2*alpha;
%     else
%         alpha=alpha0;
%     end
% elseif mu>3e-5 && iter>=5  &&  mu<3e-4
%     alpha0=0.25;
%      if alpha<=alpha0/4
%         alpha=2*alpha;
%     else
%         alpha=alpha0;
%     end
% elseif  mu<3e-6  &&  phi>3e-7
%     alpha0=0.5;
%      if alpha<=alpha0/4
%         alpha=2*alpha;
%     else
%         alpha=alpha0;
%     end
% elseif  mu<3e-6  &&  phi<3e-7 
%     alpha0=1;
%     if alpha<=alpha0/4
%         alpha=2*alpha;
%     else
%         alpha=alpha0;
%     end
% else
%     alpha=1;
% end
    %alpha=1;
    linsh=1;
    newx=x+alpha*dx;
    news=s+alpha*ds;
    newla=lambda+alpha*dla;
    nmu = mu + alpha*du;
    [newz,newt] = compz(rho,nmu,newx,news);
    nres=newz-newx;
    nphi=0.5*(nres'*nres);

    while nphi>=(1-2*sigma*alpha)*phi 
        alpha=delta*alpha;
        linsh=linsh+1;
        newx=x+alpha*dx;
        news=s+alpha*ds;
        newla=lambda+alpha*dla;
        nmu = mu + alpha*du;
        [newz,newt] = compz(rho,nmu,newx,news);
        nres=newz-newx;
        %nres=[matTvec(input_m,input_n,newla,Spt)+news-c_Spt;newz-newx;matvec(newx,Spt)-b0];
        %nres=[newz-newx;matvec(newx,Spt)-b0];
        nphi=0.5*(nres'*nres);
        
        if linsh>10
            disp('it has too much line searches')     %alpha太小！
            break 
        end
    end
    iter=iter+1;
    rho0=rho;
    rho=max(rho0,beta*norm(news,inf)/max(norm(newx),1));
    if iter>15  &&  rho>100*rho0 
        opt=0;
        disp('primal nondegenerate condition not holds ')
        break
    else
    x=newx;
    s=news;
    lambda=newla;    
    end

    if rho>rho0
        [z,t] = compz(rho,nmu,x,s);
        %res=[matTvec(input_m,input_n,lambda)+s-c;z-x;matvec(input_m,input_n,x)-b];
        %res=[z-x;matvec(input_m,input_n,x)-b0];
        res=z-x;
        phi=0.5*(res'*res);
    else
        phi=nphi;
        z=newz;
        t=newt;
    end
    mu=nmu;
    %y=mu./z;
    
    while mu > max(epsi,eta*phi) 
        mu = mu /eta;
        [z,t] = compz(rho,mu,x,s);
        %res=[matTvec(input_m,input_n,lambda)+s-c;z-x;matvec(input_m,input_n,x)-b];
        %res=[z-x;matvec(input_m,input_n,x)-b0];  
        res=z-x;
        phi=0.5*(res'*res);
    end
    output(pl,iter,norm_rh,norm_rd,phi,mu,c_Spt'*x,b0'*lambda,rho,alpha);
end

u=c_Spt-matTvec_Spt(input_m,input_n,lambda, Spt_i,Spt_j);
u=u(u<0);
pinf= norm(matvec_Spt0(x,Spt,Spt_i,Spt_j)-b)/(1+norm(b));
gap = abs(c_Spt'*x-b0'*lambda)/(1+abs(c_Spt'*x)+abs(b0'*lambda));
%dinf= norm(u)/(1+norm(x)+norm(s));
dinf=norm(u)/(1+norm(c_Spt));
opt=[opt;c_Spt'*x;b0'*lambda];
residue=[pinf;gap;dinf];

fprintf('iterations: %4d\n', iter);
fprintf('primal feasibility: %8.2e\n', pinf);
fprintf('dual feasibility: %8.2e\n', norm(matTvec_Spt(input_m,input_n,lambda, Spt_i,Spt_j)+s-c_Spt)/(1+norm(c_Spt)));
fprintf('duality gap: %8.2e\n', gap);
fprintf('dinf: %8.2e\n', dinf);
end



function header(pl)
    if (pl >= 1)
        fprintf(' ');
        fprintf('%4s    ', 'iter');
        fprintf('%8s  ', 'pr feas');
        fprintf('%8s  ', 'dl feas');
        fprintf('%8s  ', ' phi');
        fprintf('%8s  ', 'mu');
        fprintf('%8s  ', '<c,x>');
        fprintf('%8s  ', '<b,lambda>');
        fprintf('%8s  ', 'rho');
        fprintf('%8s  ', 'alpha');
    end
    if (pl >= 1)
        fprintf('\n ====    ========  ========   ========   ========  ========');
    end
    if (pl >= 2)
        fprintf('  ========  ========  ========');
    end
    if (pl >= 1) 
        fprintf('\n'); 
    end
end


function output(pl,it,xinf,sinf,rinf,mu,Pobj,Dobj,rho,alpha)
    if (pl >= 1)
        fprintf(' ');
        fprintf('%4d    ', it);
        fprintf('%8.2e  ', xinf);
        fprintf('%8.2e  ', sinf);
        fprintf('%8.2e  ', rinf);
        fprintf('%8.2e  ', mu);
        fprintf('%8.2e  ', Pobj);
        fprintf('%8.2e  ', Dobj);
        fprintf('%8.2e  ', rho);
        fprintf('%8.2e  ',alpha);
    end

    if (pl >= 1) 
        fprintf('\n');
    end
end


 function result = matvec_Spt0(x,Spt,Spt_i,Spt_j)
  X = sparse(Spt_i, Spt_j, x, size(Spt,1), size(Spt,2));
  x1=sum(X,1)';
  x2=sum(X,2);
  result=[x1;x2];
 end

  function result = matvec_Spt(x,Spt,Spt_i,Spt_j)
  X = sparse(Spt_i, Spt_j, x, size(Spt,1), size(Spt,2));
  x1=sum(X,1)';
  x2=sum(X,2);
  x2=x2(1:end-1);
  result=[x1;x2];
  end

  
function result = matTvec_Spt(input_m, input_n, lambda, Spt_i,Spt_j)      
  lambda=[lambda;0];
  lambda1=lambda(1:input_m);
  lambda2=lambda(input_n+1:end);
  result=lambda1(Spt_j)+lambda2(Spt_i);
end 

function result = MAIN_FUNC(D,V2,p,x)
    result = (V2(p)).*x - D*(D'*x);
end


function result = MAIN_FUNC2(DD,V1,V2,x)
    result = (V2).*x - DD*((DD'*x)./(V1'));
end


function result = PRE_FUN2(V2,x)
    %result = x./(sqrt(V2));
     result = x./(V2);
end

function [z,t] = compz(rho,mu,x,s)
	v = s-rho.*x;
    t = sqrt(v.^2+4*mu*rho);
	z = (t - v)./(2*rho);
end

function [x,s,lambda]= Warmstart_Spt(m,n,b,c_Spt,x,lambda,s,mu,rho,R,p,Spt,Spt_i,Spt_j) 
	    maxIter =10; 
        b0=b(1:end-1);
        V1 = sum(Spt,1)'+1e-12;     
        St=Spt(1:end-1, :);
        z = compz(rho,mu,x,s);
	for iter = 1:maxIter
		rhs = matvec_Spt(c_Spt-rho*z-s,Spt,Spt_i,Spt_j)+rho*b0;
        rhs1 = rhs(1:m);
        rhs2 = rhs(m+1:m+n-1);
        rhsSchur0 = rhs2 - St*(rhs1./V1);
        la2=zeros(size(rhsSchur0));
        la2(p) = R\(R'\(rhsSchur0(p)));
        la1 = (rhs1 - St'*la2) ./ V1;
        lambda= [la1; la2];
		x = z + (matTvec_Spt(m, n, lambda, Spt_i,Spt_j)+s-c_Spt)/rho;
		z = compz(rho,mu,x,s);
        if norm(z-x) < 1e-2
            return
        end
		s = s + rho*(z-x);
	end
end