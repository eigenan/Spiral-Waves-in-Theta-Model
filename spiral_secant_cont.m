clear all;                          close all;
tic
%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_init=0.0;                        mu_end=2.0;
r_min=0.1;                          r_max=70;
p_min=0;                            p_max=2*pi;
r_grid_size=600;                    p_grid_size=60;
dr=(r_max-r_min)./r_grid_size;      dp=(p_max-p_min)./p_grid_size;
r=(r_min:dr:r_max)';                Nr=length(r);
p=(p_min:dp:p_max)';                Np=length(p);
n=Nr*Np;                            pars=[n Nr Np];
g=@(u,mu) mu.*cos(u)+1;             dg=@(u,mu) -mu.*sin(u);
dg_mu=@(u,mu) cos(u);               tol=1e-6;
ds=5;                               increment_init=1e-1;
mu_tol=1e-2;                        overshoot=0.7;
few_newton=1.3;                     newton_fail=0.7;
maxit_newton=6;                     maxit_gmres=6;
mu_list=[];                         wavenumber=[];
wnp=floor(p_grid_size/2);           wn_trim=floor(0.8*r_grid_size);
omega_list=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
er=ones(Nr,1);                  ep=ones(Np,1);
% 1/r % 1/r^2
Lr=kron(spdiags(1./r,0,Nr,Nr),speye(Np,Np));
Lrr=kron(spdiags(1./r.^2,0,Nr,Nr),speye(Np,Np));
% Neumann & periodic boundary condition
neu_r=sparse([1 Nr],[2 Nr-1],[-0.5 0.5],Nr,Nr);
neu_rr=sparse([1 Nr],[2 Nr-1],[1 1],Nr,Nr);
neu_p=sparse([1 Np],[Np 1],[-0.5 0.5],Np,Np);
neu_pp=sparse([1 Np],[Np 1],[1 1],Np,Np);
% dr
DDr=(spdiags([-er./2 zeros(Nr,1) er./2],-1:1,Nr,Nr)+neu_r)/dr;
Dr=kron(DDr,speye(Np));
% dp
DDp=(spdiags([-ep./2 zeros(Np,1) ep./2],-1:1,Np,Np)+neu_p)/dp;
Dp=kron(speye(Nr),DDp);

% neu_r=sparse([1 Nr],[1 Nr],[-1 1],Nr,Nr);
% neu_rr=sparse([1 Nr],[1 Nr],[1 1],Nr,Nr);
% neu_p=sparse([1 Np],[Np 1],[-1 1],Np,Np);
% neu_pp=sparse([1 Np],[Np 1],[1 1],Np,Np);
% % dr
% DDr=(spdiags([-er zeros(Nr,1) er],-1:1,Nr,Nr)+neu_r)/dr/2;
% Dr=kron(DDr,speye(Np));
% % dp
% DDp=(spdiags([-ep zeros(Np,1) ep],-1:1,Np,Np)+neu_p)/dp/2;
% Dp=kron(speye(Nr),DDp);

% laplacian (second order)
lapr=(spdiags([er -2*er er],-1:1,Nr,Nr)+neu_rr)/dr^2;
lapp=(spdiags([ep -2*ep ep],-1:1,Np,Np)+neu_pp)/dp^2;
% full equation f(w)
A=@(w) -w(end).*Dp + (kron(lapr,speye(Np)) + Lr * Dr + Lrr * kron(speye(Nr),lapp));
f=@(w,mu) [A(w)*w(1:end-1) + g(w(1:end-1)+reshape(p*ones(1,Nr),n,1),mu) - w(end).*ones(n,1); sum(w(1:end-1))];
f_mu=@(nu,nu_fixed,ek) [f(nu(1:end-1,1),nu(end,1)); dot((nu-nu_fixed),ek)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial Solutions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu0=sparse([n+1 n+2],[1 1],[1+mu_init mu_init],n+2,1);
nu1=nu0;
error=1;
counter=1;
while error > tol && counter < maxit_newton
    res = f(nu0(1:end-1,1),nu1(end,1));
    nu1(1:end-1,1) = nu0(1:end-1,1) - df(nu0(1:end-1,1),dg,A,Dp,p,nu0(end,1),pars)\res;
    error=norm(res);
    nu0(1:end-1,1)=nu1(1:end-1,1);
    counter=counter+1;
end
% plot solution
figure(1)
v=reshape(nu1(1:end-2,1),Np,Nr);
u=mod(v+p*ones(1,Nr),2*pi);
imagesc(r,p,u)
colorbar
title(['mu=' num2str(nu1(end,1)) ', frequency=' num2str(nu1(end-1,1)) ', homog osc freq =' num2str(sqrt(1-nu1(end,1).^2))])
xlabel("$r$","interpreter","latex","FontSize",18)
ylabel("$\varphi$","interpreter","latex","FontSize",18)
drawnow

figure(2)
[theta,rho]=ndgrid(p,r);
[x,y]=pol2cart(theta,rho);
surf(x,y,u,'EdgeColor','none')
% title(['u mod 2\pi with mu=' num2str(nu(end,1))])
title("$u$ mod $2\pi$","interpreter","latex","FontSize",24)
xlabel("$x$","interpreter","latex","FontSize",24)
ylabel("$y$","interpreter","latex","FontSize",24)
pbaspect([1 1 1])
view(2)

nu_buffer=nu0;
nu1=nu0+sparse(n+2,1,increment_init,n+2,1);

if nu1(end,1) < mu_end
    error=1;
    counter=1;
    while error > tol && counter < maxit_newton
        res = f(nu0(1:end-1,1),nu1(end,1));
        nu1(1:end-1,1) = nu0(1:end-1,1) - df(nu0(1:end-1,1),dg,A,Dp,p,nu0(end,1),pars)\res;
        error=norm(res);
        nu0(1:end-1,1)=nu1(1:end-1,1);
        counter=counter+1;
    end
    % plot solution
    figure(1)
    v=reshape(nu1(1:end-2,1),Np,Nr);
    u=mod(v+p*ones(1,Nr),2*pi);
    imagesc(r,p,u)
    colorbar
    title(['mu=' num2str(nu1(end,1)) ', frequency=' num2str(nu1(end-1,1)) ', homog osc freq =' num2str(sqrt(1-nu1(end,1).^2))])
    xlabel("$r$","interpreter","latex","FontSize",18)
    ylabel("$\varphi$","interpreter","latex","FontSize",18)
    drawnow

    figure(2)
    [theta,rho]=ndgrid(p,r);
    [x,y]=pol2cart(theta,rho);
    surf(x,y,u,'EdgeColor','none')
    % title(['u mod 2\pi with mu=' num2str(nu(end,1))])
    title("$u$ mod $2\pi$","interpreter","latex","FontSize",24)
    xlabel("$x$","interpreter","latex","FontSize",24)
    ylabel("$y$","interpreter","latex","FontSize",24)
    pbaspect([1 1 1])
    view(2)
end

nu0=nu_buffer;
nu=nu1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Secant Continuation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
halt=0;
counter=1;
while nu(end,1)<mu_end

    halt=0;
    nu_buffer=nu;
    ek = (nu1 - nu0)./norm(nu1 - nu0);
    nu_fixed = nu1 + ds.*ek;
    nu=nu_fixed;

    error=1;
    counter=1;
    while error > tol && counter < maxit_newton

        %%%%%%%%%%%%%%%%%%%%%% x=A\b %%%%%%%%%%%%%%%%%%%%%%%%
        res=f_mu(nu,nu_fixed,ek);
        nu_new = nu - jacF(nu,ek,dg,dg_mu,A,Dp,p,pars)\res;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%% GMRES %%%%%%%%%%%%%%%%%%%%%%%%
        % jacf=df(w1,dg,A,Dp,p,mu,pars);
        % [L,U] = ilu(jacf,struct('type','ilutp','droptol',1e-6));
        % lin_sol = gmres(jacf,f(w1,mu),3,tol,maxit_gmres,L,U);
        % w = w1 - lin_sol;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        error=norm(res);
        nu=nu_new;
        counter=counter+1;
    end

    % plot solution
    v=reshape(nu(1:end-2),Np,Nr);
    u=mod(v+p*ones(1,Nr),2*pi);
    figure(1)
    imagesc(r,p,u)
    colorbar
    % title(['mu=' num2str(nu(end,1)) ', frequency=' num2str(nu(end-1,1)) ', homog osc freq =' num2str(sqrt(1-nu(end,1).^2))])
    % title(['u mod 2\pi with mu=' num2str(nu(end,1))])
    title("$u$ mod $2\pi$","interpreter","latex","FontSize",24)
    xlabel("$r$","interpreter","latex","FontSize",24)
    ylabel("$\varphi$","interpreter","latex","FontSize",24)
    drawnow

    figure(2)
    [theta,rho]=ndgrid(p,r);
    [x,y]=pol2cart(theta,rho);
    surf(x,y,u,'EdgeColor','none')
    % title(['u mod 2\pi with mu=' num2str(nu(end,1))])
    title("$u$ mod $2\pi$","interpreter","latex","FontSize",24)
    xlabel("$x$","interpreter","latex","FontSize",24)
    ylabel("$y$","interpreter","latex","FontSize",24)
    pbaspect([1 1 1])
    view(2)

    disp("mu="+num2str(nu(end,1)))
    disp("omega="+num2str(nu(end-1,1)))

    if nu(end,1) > mu_end + mu_tol
        ds = overshoot.*ds;
        disp("Overshoot. Decrease step size to ds="+ds)
        nu=nu_buffer;
    else
        if (counter >= maxit_newton && error >= tol)
            halt=1;
            ds = newton_fail.*ds;
            disp("Newton failed. Decrease step size to ds="+ds)
        else
            halt=0;
            nu0=nu1;
            nu1=nu;
            %%%%%%%%%%%%%%%%%%%%%%%% wavenumber %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(mu_list) | nu(end,1) > mu_list(end) + 0.01
                pos = wave_front(trimdata(u(wnp,:),wn_trim,Side="both"));
                dist=0;
                if length(pos) >= 2
                    for j=1:length(pos)-1
                        dist = dist + abs(pos(j+1) - pos(j));
                    end
                    dist = dist/(length(pos)-1);
                end
                if dist ~= 0
                    mu_list(end+1) = nu(end,1);
                    wavenumber(end+1) = 2*pi./(dr.*dist);
                    omega_list(end+1) = nu(end-1,1);
                    disp("mu="+mu_list(end)+" wavenumber="+wavenumber(end)+" omega="+omega_list(end))
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    
        if counter < maxit_newton && error < tol
            ds = few_newton.*ds;
            disp("Too few Newton steps. Increase step size to ds="+ds)
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot wavenumber-mu graph
figure(3)
plot(mu_list,wavenumber,'b.-','LineWidth',2)
title("$\mu$ vs the wavenumber","interpreter","latex","FontSize",24)
xlabel("$\mu$","interpreter","latex","FontSize",24)
ylabel("wavenumber","interpreter","latex","FontSize",24)

% plot omega-mu graph
figure(4)
plot(mu_list,omega_list,'r.-','LineWidth',2)
title("$\mu$ vs $\omega$","interpreter","latex","FontSize",24)
xlabel("$\mu$","interpreter","latex","FontSize",24)
ylabel("$\omega$","interpreter","latex","FontSize",24)

disp('Done!')
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jacF=jacF(nu1,ek,dg,dg_mu,A,Dp,p,pars)
    % Jacobian of F(w,mu), boardered matrix with an extra argument mu
    n=pars(1);
    Nr=pars(2);
    N=length(nu1);
    mu=nu1(end,1);
    jacF=sparse(N,N);
    jacF(1:end-1,1:end-1) = df(nu1(1:end-1,1),dg,A,Dp,p,mu,pars);
    jacF(1:end-2,end) = dg_mu(nu1(1:end-2,1)+reshape(p*ones(1,Nr),n,1),nu1(end,1));
    jacF(end,1:end)=ek';
end

function df=df(w,dg,A,Dp,p,mu,pars)
    % Jacobian Df(w)
    n=pars(1);
    Nr=pars(2);
    N=length(w);
    df=sparse(N,N);
    dg_block=sparse(1:n,1:n,dg(w(1:end-1)+reshape(p*ones(1,Nr),n,1),mu),n,n);
    df(1:end-1,1:end-1) = A(w) + dg_block;
    df(1:end-1,end) = -Dp*w(1:end-1) - ones(n,1);
    df(end,1:end-1) = ones(1,n);
end

function pos=wave_front(v)
    % input: v a row vector (a ray in radial direction)
    % output: a vector consisting of indices of elements at a jump of almost 2*pi
    pos = [];
    for j=1:length(v)-1
        if abs(v(j+1)-v(j)) > 1.9*pi
            pos(end+1) = j;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%