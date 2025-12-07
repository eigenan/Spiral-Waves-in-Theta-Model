clear all;                          close all;
vid = VideoWriter('spiral.mp4','MPEG-4');
open(vid)
%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_over=0;                       mu=0.99;
dt=1e-1;                            T=100;   
D2=1;                               D4=0; % for unstable spirals: set D2 as a small negative number and turn on D4
g=@(u,mu) mu.*cos(u)+1;             dg=@(u,mu) -mu.*sin(u);
r_min=0.1;                          r_max=80;
p_min=0;                            p_max=2*pi;
r_grid_size=1000;                    p_grid_size=120;
dr=(r_max-r_min)./r_grid_size;      dp=(p_max-p_min)./p_grid_size;
r=(r_min:dr:r_max)';                Nr=length(r);
p=(p_min:dp:p_max)';                Np=length(p);
n=Nr*Np;                            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

er=ones(Nr,1);                  ep=ones(Np,1);
% 1/r % 1/r^2
Lr=kron(spdiags(1./r,0,Nr,Nr),speye(Np,Np));
Lrr=kron(spdiags(1./r.^2,0,Nr,Nr),speye(Np,Np));
% Neumann & periodic boundary condition
bc_r=sparse([1 Nr],[2 Nr-1],[-0.5 0.5],Nr,Nr);
bc_rr=sparse([1 Nr],[2 Nr-1],[1 1],Nr,Nr);
bc_p=sparse([1 Np],[Np 1],[-0.5 0.5],Np,Np);
bc_pp=sparse([1 Np],[Np 1],[1 1],Np,Np);
% dr
DDr=(spdiags([-er./2 zeros(Nr,1) er./2],-1:1,Nr,Nr)+bc_r)/dr;
Dr=kron(DDr,speye(Np));
% dp
DDp=(spdiags([-ep./2 zeros(Np,1) ep./2],-1:1,Np,Np)+bc_p)/dp;
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
lapr=(spdiags([er -2*er er],-1:1,Nr,Nr)+bc_rr)/dr^2;
lapp=(spdiags([ep -2*ep ep],-1:1,Np,Np)+bc_pp)/dp^2;

lap=(kron(lapr,speye(Np)) + Lr * Dr + Lrr * kron(speye(Nr),lapp));
A=(speye(n,n)-dt*(-D4*lap*lap+D2*lap));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=0;
if start_over==1
    v0=zeros(n,1);
else
    load("dns_data.mat");
    v0=v;
end


v=reshape(v0,Np,Nr);
u=mod(v+p*ones(1,Nr),2*pi);
% u=v+p*ones(1,Nr);
[theta,rho]=ndgrid(p,r);
[x,y]=pol2cart(theta,rho);
surf(x,y,u,'EdgeColor','none')
view(2)
% imagesc(r,p,u)
colorbar
title("t="+t)
xlabel("$r$","interpreter","latex","FontSize",18)
ylabel("$\varphi$","interpreter","latex","FontSize",18)
axis square
drawnow

while t<T

    v=A\(v0+dt*g(v0+reshape(p*ones(1,Nr),n,1),mu));
    % if max(v) > 5
    %     v=0.5.*v;
    % end
    v0=v;
    t=t+dt;
    save("dns_data.mat","v");

    v=reshape(v,Np,Nr);
    u=mod(v+p*ones(1,Nr),2*pi);
    % u=v+p*ones(1,Nr);
    [theta,rho]=ndgrid(p,r);
    [x,y]=pol2cart(theta,rho);
    surf(x,y,u,'EdgeColor','none')
    view(2)
    % imagesc(r,p,u)
    colorbar
    title("t="+t)
    xlabel("$r$","interpreter","latex","FontSize",18)
    ylabel("$\varphi$","interpreter","latex","FontSize",18)
    axis square
    drawnow
    frame = getframe(gcf);
    writeVideo(vid,frame);

end

close(vid)

disp('Done!')