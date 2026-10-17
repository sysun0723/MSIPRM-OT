function [x,s,lambda,z,opt,iter] = iprm4ot0(input_m,input_n,Su,De,printlevel)
persistent use
tol=1e-6;
maxit=100;
pl = printlevel;

C=construct_C(sqrt(input_m),sqrt(input_n));
c=C(:);
Su = Su/sum(Su);
De = De/sum(De); 
b = [Su; De];
b=b/norm(b);
b0=b(1:end-1);


[x,lambda,s] = warmstart(input_m,input_n,c,b);

opt = -1;
rho=10;        %balance parameter, current not update
gamma0=1e-3;
epsi=tol;    %stopping parameter
delta=0.5;    %stepsize parameter
sigma=0.01;   %line search parameter
%maxitpcg = 500;
%mu=max(abs(x).*abs(s));
%mu = max(mu,1);
mu =1;
eta = 5;
[z,t] = compz(rho,mu,s,x);
%y=mu./z;
res=[matTvec0(input_m,input_n,lambda)+s-c;z-x;matvec0(input_m,input_n,x)-b0];
% z=zeros(n,1);newz=zeros(n,1);
% y=zeros(n,1);newy=zeros(n,1);
%  z(pos_vars)=0.5*(sqrt((s(pos_vars)-rho*x(pos_vars)).^2+4*mu*rho)-(s(pos_vars)-x(pos_vars)*rho))/rho;
%  y(pos_vars)=0.5*(sqrt((s(pos_vars)-rho*x(pos_vars)).^2+4*mu*rho)+(s(pos_vars)-x(pos_vars)*rho))/rho;
% res=[matTvec(input_m,input_n,lambda)+s-c;z(pos_vars)-x(pos_vars);matvec(input_m,input_n,x)-b];
phi=0.5*(res'*res);
norm_b=norm(b0);
norm_c=norm(c);
while mu > eta*max(epsi/eta,phi)
    mu = mu /eta;
    [z,t] = compz(rho,mu,s,x);
    %res=[matTvec0(input_m,input_n,lambda)+s-c;z-x;matvec0(input_m,input_n,x)-b0]; 
    %y(pos_vars)=0.5*(sqrt((s(pos_vars)-rho*x(pos_vars)).^2+4*mu*rho)+(s(pos_vars)-x(pos_vars)*rho))/rho;
    %z(pos_vars) = (2*mu)./(sqrt((s(pos_vars)-rho*x(pos_vars)).^2+4*mu*rho)+(s(pos_vars)-x(pos_vars)*rho));
    %res=[matTvec(input_m,input_n,lambda)+s-c;z-x;matvec(input_m,input_n,x)-b]; 
    %res=[z-x;matvec(input_m,input_n,x)-b]; 
    res=z-x; 
    phi=0.5*(res'*res);
end
gamma = min(gamma0,mu/phi);
iter=0;
header(pl); 
while iter<maxit 
    rh = matvec0(input_m,input_n,x)-b0;
    rd =c-matTvec0(input_m,input_n,lambda)-s;
    norm_res=[abs(c'*x-b0'*lambda)/(1+abs(c'*x)+abs(b0'*lambda));norm(rh)/(1+norm_b);norm(rd)/(1+norm_c);mu];
    % if (nnz(isnan(rh)) > 0 || nnz(isinf(rh)) > 0 || nnz(isnan(rd)) > 0 || nnz(isinf(rd)) > 0)
    %     opt = -1;
    %     return;
    % end
    if isempty(find(norm_res >= epsi,1)) 
        disp('problem solved');
%        [norm(rh)/(1+norm_b);norm(rd)/(1+norm_c);abs(c'*x-b0'*lambda)/(1+abs(c'*x))];
        opt=1;
        break;
    end    
    
    %The search direction
%     norm_rb=norm(rh);
%     norm_rc=norm(rd);
% if 1e6*norm_rb<mu || 1e6*norm_rc< mu
% du =-mu +0.5*gamma*phi;
% elseif 1e4*norm_rb<mu || 1e4*norm_rc<mu
% du =-mu +0.5*gamma*phi;
% elseif 1e3*norm_rb<mu || 1e3*norm_rc <mu
% du=-mu +0.75*gamma*phi;
% elseif 1e2*norm_rb<mu || 1e2*norm_rc<mu
% du =-mu +0.75*gamma*phi;
% else
% du=-mu+0.995*gamma*phi;
% end    
if mu>3e-4  &&  iter<6
    gamma0=5e-3;
    du=-mu+gamma*phi;
elseif mu>3e-4 && iter>=5 
    gamma0=5e-3;
    du=-mu+gamma*phi;
elseif mu>3e-5 && iter>=5  &&  mu<3e-4
    gamma0=5e-3;
    du=-mu+gamma*phi;
elseif  mu<3e-6  &&  phi>1e-8
    gamma0=1e-3;
    du=-mu+gamma*phi;
elseif  mu<3e-6  &&  phi<1e-8  && phi>1e-9
    gamma0=1e-3;
    du=-mu+0.5*gamma*phi;
elseif  mu<3e-6  &&  mu>1e-8   &&  phi<1e-9
    gamma0=3e-3;
    du=-mu+0.8*gamma*phi;
else
    gamma0=2e-3;
    du=-mu+gamma*phi;
end

    rf = du+(t.*(z-x));
    rhs = mu*b0-matvec0(input_m,input_n,mu*x+z.*(-z.*rd+rf));
    D = reshape(z.*z,[input_m,input_n]);
    V1 = sum(D,1)';   
    V2 = sum(D,2);
    V2 = V2(1:end-1);  
    DD=D(1:end-1, :);
        rhs1 = rhs(1:input_n);                  
        rhs2 = rhs(input_n+1:end);                
        rhsSchur = rhs2 - DD*(rhs1./V1);         
        DD0 = DD'./sqrt(V1);                
        main_func3 = diag(V2) - DD0'*DD0;
        if isempty(use)
        use= 1;  
        end  
        if use 
           try 
            R = chol(main_func3);
            dla2 = R\(R'\rhsSchur);
           catch 
            use= 0; 
            [L, D, P] = ldl(main_func3, 'vector');
            warning('off','all');
            dla2 = L'\(D\(L\rhsSchur(P)));
           end 
        else 
        [L, D, P] = ldl(main_func3, 'vector');
        warning('off','all');
        dla2 = L'\(D\(L\rhsSchur(P)));
        end 
        dla1 = (rhs1 - DD'*dla2)./V1;      
        dla = [dla1; dla2];                     
    ds = rd-matTvec0(input_m,input_n,dla);
    %dx = (re-z.*ds)./(rho*y);
    dx = (z.*(rf-z.*ds))./mu;
    %Doing line search

if mu>3e-4  &&  iter<6
    alpha=1;
elseif mu>3e-4 && iter>=5 
    alpha0=0.0625;
    if alpha<=alpha0/8
        alpha=2*alpha;
    else
        alpha=alpha0;
    end
elseif mu>3e-5 && iter>=5  &&  mu<3e-4
    alpha0=0.125;
    if alpha<=alpha0/8
        alpha=2*alpha;
    else
        alpha=alpha0;
    end
elseif  mu<3e-5  &&  mu>3e-6
    alpha0=0.25;
    if alpha<=alpha0/8
        alpha=2*alpha;
    else
        alpha=alpha0;
    end
elseif  mu<3e-6  && mu>1e-6
    alpha0=0.5;
    if alpha<=alpha0/8
        alpha=2*alpha;
    else
        alpha=alpha0;
    end
else
    alpha=1;
end
    linsh=1;
    newx=x+alpha*dx;
    news=s+alpha*ds;
    newla=lambda+alpha*dla;
    nmu = mu + alpha*du;
    [newz,newt] = compz(rho,nmu,news,newx);
    nres=newz-newx;
    %nres=[matTvec0(input_m,input_n,newla)+news-c;newz-newx;matvec0(input_m,input_n,newx)-b0];
    %nres=[newz-newx;matvec0(input_m,input_n,newx)-b];
    nphi=0.5*(nres'*nres);

    while nphi>=(1-2*sigma*alpha)*phi 
        alpha=delta*alpha;
        linsh=linsh+1;
        newx=x+alpha*dx;
        news=s+alpha*ds;
        newla=lambda+alpha*dla;
        nmu = mu + alpha*du;
        [newz,newt] = compz(rho,nmu,news,newx);
        %newz(pos_vars)=0.5*(sqrt((news(pos_vars)-rho*newx(pos_vars)).^2+4*nmu*rho)-(news(pos_vars)-newx(pos_vars)*rho))/rho;
        %newy(pos_vars)=0.5*(sqrt((news(pos_vars)-rho*newx(pos_vars)).^2+4*nmu*rho)+(news(pos_vars)-newx(pos_vars)*rho))/rho;
%         if (nnz(isnan(newz)) > 0 || nnz(isinf(newz)) > 0 )
%             newz = (2*nmu)./(sqrt((news-rho*newx).^2+4*nmu*rho)+(news-newx*rho));
%         end
        %nres=[newz-newx;matvec0(input_m,input_n,newx)-b0];
        nres=newz-newx;
        %nres=[matTvec0(input_m,input_n,newla)+news-c;newz-newx;matvec0(input_m,input_n,newx)-b0];
        nphi=0.5*(nres'*nres);
        
        if linsh>20
            disp('it has too much line searches')     %alpha太小！
            break 
        end
    end
    
    if  nphi>1 && iter>15
       break
    else
    x=newx;
    s=news;
    lambda=newla;
    end
    iter=iter+1;
    rho0=rho;
    %rho=max(rho0,0.1*norm(s,inf)/max(norm(x),1));
    rho=max(rho0,0.1*norm(s,inf)/max(norm(x),1));
    if rho>rho0
        [z,t] = compz(rho,nmu,s,x);
        res=z-x;
        phi=0.5*(res'*res);
    else
        %y=newy;
        %res=nres;
        phi=nphi;
        z=newz;
        t=newt;
    end
    mu=nmu;
    %y=mu./z;

    %% 如果注释掉这块算的更快精度更高但目标值不太对
    while mu > eta*max(epsi/eta,phi)
        mu = mu /eta;
        [z,t] = compz(rho,mu,s,x);
        %z(pos_vars)=0.5*(sqrt((s(pos_vars)-rho*x(pos_vars)).^2+4*mu*rho)-(s(pos_vars)-x(pos_vars)*rho))/rho;
        %y(pos_vars)=0.5*(sqrt((s(pos_vars)-rho*x(pos_vars)).^2+4*mu*rho)+(s(pos_vars)-x(pos_vars)*rho))/rho;
%         if (nnz(isnan(z)) > 0 || nnz(isinf(z)) > 0 )
%              z = (2*mu)./(sqrt((s-rho*x).^2+4*mu*rho)+(s-x*rho));
%         end
        %res=[matTvec(input_m,input_n,lambda)+s-c;z-x;matvec(input_m,input_n,x)-b];
        %res=[z-x;matvec0(input_m,input_n,x)-b0];  
        res=z-x;
        phi=0.5*(res'*res);
    end
    gamma = min(gamma0,mu/phi);
    output(pl,iter,norm(rh),norm(rd),phi,mu,c'*x,b0'*lambda,rho,alpha);
end

u=c-matTvec0(input_m,input_n,lambda);
u=u(u<0);
uu=x-s;
uu(uu<0)=0;
uu=x-uu;  

fprintf('iterations: %4d\n', iter);
fprintf('primal feasibility: %8.2e\n', norm(matvec0(input_m,input_n,x)-b0)/(1+norm(b0)));
fprintf('dual feasibility: %8.2e\n', norm(matTvec0(input_m,input_n,lambda)+s-c)/(1+norm(c)));
fprintf('duality gap: %8.2e\n', abs(c'*x-b0'*lambda)/(1+abs(c'*x)+abs(b0'*lambda)));
% fprintf('complementarity: %8.2e\n', abs(full(dot(x,s))/N));  
% fprintf('relative z-x error: %8.2e\n', norm(z-x)/(1+norm(x)+norm(z)));
% fprintf('relative optimal error: %8.2e\n', abs(c'*(z-x))/(1+abs(c'*z)+abs(c'*x)));
 fprintf('dinf: %8.2e\n', norm(u)/(1+norm(c)));
% fprintf('KKT residue: %8.2e\n', norm(uu)/(1+norm(x)+norm(s)));
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

 function result = matvec0(input_m,input_n,x)
 P=reshape(x,[input_m,input_n]);
 x1=sum(P,1)';
 x2=sum(P,2);
 xx=[x1;x2];
 result=xx(1:end-1);
 end
 
 function result = matTvec0(input_m,input_n,lambda)
  lambda=[lambda;0];
  lambda1=lambda(1:input_n);
  lambda2=lambda(input_n+1:input_n+input_m)';
  Res=(lambda1+lambda2)';
  result=Res(:);
 end

function [z,t] = compz(rho,mu,s,x)
	v = s-rho.*x;
    t = sqrt(v.^2+4*mu*rho);
	z = (t - v)./(2*rho);
%     n = size(s,1);
%     ind1 = find(v < 0);
%     ind2 = setdiff((1:n)',ind1);
%     z = zeros(n,1);
%     z(ind1) = (t(ind1) - v(ind1))./(2*rho);
%     z(ind2) = (2*mu)./(t(ind2)+v(ind2));
 %   z = (2*mu)./(t + v);
end

function [x,lambda,s] = warmstart(input_m,input_n,c,b)
    b=b(1:end-1);
	rho = 10;
	mu = 0.01;
	maxIter =50;  
    b1 = b(1:input_n);
    b2 = b(input_n+1:end);
    bb=b2-sum(b1,1)./input_m;

    x2=bb./input_n+sum(bb,1)/input_n;
    x1=b1./input_m-sum(x2,1)/input_m;
    xx=[x1;x2];
    x= matTvec0(input_m,input_n,xx);
    
    c0=matvec0(input_m,input_n,c);
    c1 = c0(1:input_n);
    c2 = c0(input_n+1:end);
    cc=c2-sum(c1,1)./input_m;

    lambda2=cc./input_n+(sum(cc,1))/(input_n);
    lambda1=c1./input_m-sum(lambda2,1)/input_m;
    lambda=[lambda1;lambda2];
    s = c - matTvec0(input_m,input_n,lambda);

	z = compz(rho,mu,s,x);
	for iter = 1:maxIter
		rhs = matvec0(input_m,input_n,c-rho*z-s)+rho*b;
        rhs1 = rhs(1:input_n);
        sum1=sum(rhs1,1);
        rhs2 = rhs(input_n+1:end);
        rhss=rhs2-sum1/input_m;
        lambda2=rhss./input_n+(sum(rhss,1))/input_n;
        lambda1=rhs1./input_m-sum(lambda2,1)/input_m;
        lambda=[lambda1;lambda2];
		x = z + (matTvec0(input_m,input_n,lambda)+s-c)/rho;
		z = compz(rho,mu,s,x);
        if norm(z-x) < 1e-3
            return
        end
		s = s + rho*(z-x);
	end
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