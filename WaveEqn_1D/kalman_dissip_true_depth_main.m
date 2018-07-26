%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PROGRAM DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program: 
% 1. solves the following 1d coupled ODE: 
% dh/dt = a*dq/dx 
% dq/dt = b*dh/dx 
% where h and q are vectors in this case 
% h is the tsunami wave height [m] 
% q is the horizontal flux [10^3*m^2/s] 
% 2. simulates the data according to the coupled ODE to generate synthetic data 
% (RUNGE-KUTTA INTEGRATION)
% 3. uses Kalman filter operator to predict the solution
% remark: the model is a vector in space. For each point in space (and at each time step)
% we try to solve for two unknowns (i.e., wave flux and height), but we assume that we can
% only record one component (the tsunami wave height)
% For a given time, the model is stored as follows:
% m=[q(t,x1)...q(t,xN)h(t,x1)...h(t,xN)]'
% WITH ARTIFICIAL DISSIPATION, NORMAL (UNSTAGGERED) GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DOMAIN PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
N=500; % nb of intervals
nm1=N+1;
nm=2*N+2; % nb of model points (N+1 for height and N+1 for flux)
L=550; % length of domain [km]
xmin=20; % minimum x [km] (don't set it to zero because shallow wave equation approximation breaks down)
irreg=load('earthsurface_x_y-2.mat'); % imports water bottom 
x_irreg=irreg.x_es; % irregularly sampled x-position [km]
y_irreg=irreg.y_es; % irregularly sampled y-position [km]
%get rid of overlapping locations in the x-position
y_irreg(diff(x_irreg)==0)=[];
x_irreg(diff(x_irreg)==0)=[];

x_reg=linspace(xmin,L,N+1)'; % regularly sampled x-position [km]
y_reg=-interp1(x_irreg,y_irreg,x_reg); % sea floor depth [km]
dx=(L-xmin)/N; % grid spacing [km]
x=x_reg; % easting [km]
M.g=10e-3; % [km/s^2]
M.H=y_reg; % water depth [km]
M.c=sqrt(M.g*M.H); % wave velocity [km/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT SEAFLOOR TOPOLOGY %%%%%%%%%%%%%%%%%%%%%%%
% plot seafloor depth against x
fig=1;
figure(fig)
plot(x_reg,y_reg,'LineWidth',5);
set(gca,'YDir','Reverse')
%title('Seafloor topology','Fontsize',50)
xlabel('Easting [km]','FontSize',15)
ylabel('Seafloor depth [km]','FontSize',15) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODELING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
tmax=2*L/min(M.c); % maximum recording time [s]
tmax=tmax/5;
order=6; % RK order
addpath SBPoperators
[D,Hinv,~,A,~,~,~] = SBPoperatorsRussian(N,dx,order); % SBP dissip operator
CFL=0.5;
t0=0; % initial time
dt=CFL*dx/max(M.c); % time sampling [s]
nt=ceil(tmax/dt); % number of time steps
cdiss=0.25; % artificial dissipation term
SAT=[M.c(1) M.c(end)]*Hinv; % boundary conditions
[RK,alpha]=rk_init(dt); % initialize RK coefficients
time=[t0:dt:t0+dt*(nt-1)]; % time vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RECEIVERS POSITION %%%%%%%%%%%%%%%%%%%%%%%%%
nr=50; % number of OBN
oxr=50; % x-position of first OBN (going from West -> East)  [km] 
or=ceil((oxr-xmin)/dx+1); % find grid point that corresponds to the origin (nearest neighbors)
oxr=xmin+(or-1)*dx; % find the x-coordinate of the left most receiver after nearest neighbors interpolation
dr=1; % OBN spacing [nb of grid point]
drx=dr*dx; % OBN spacing [km]
nd=nr; % number of data point (= nb of OBN) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RECEIVERS POSITION FOR QC %%%%%%%%%%%%%%%%%%
x_receivers=[oxr:drx:oxr+drx*(nr-1)];
x_receivers_grid=[or:dr:or+(nr-1)*dr];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% INITIAL CONDITIONS - LOAD GABE'S DATA %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load initial wave height at the moment of rupture
% on irregular grid
height_irreg=load('tsunami.mat');  

% load water bottom 1D profile
% on irregular grid
irreg=load('earthsurface_x_y-2.mat'); 

% define x position on irregular grid [km]
x_irreg=irreg.x_es; 

% load tsunami initial wave height 
% on irregular grid
height_irreg_load=height_irreg.tsunami;
height_irreg_load(diff(x_irreg)==0)=[];
x_irreg(diff(x_irreg)==0)=[];

% interpolate tsunami initial wave height
% on regular grid
height_reg=interp1(x_irreg,height_irreg_load,x_reg); % sea floor depth [km]
xmax=xmin+(nm1-1)*dx;

% compute weighting to remove the body waves in Gabe's simulation
for ix=1:nm1
    easting=xmin+dx*(ix-1);
    if (easting > 300)
        weight(ix,1)=sin(3.14/2*(xmax-easting)/(xmax-300))^10; 
    else
        weight(ix,1)=1;
    end
end

% apply weighting to initial wave height
height_reg_weight=height_reg.*weight;
m_true=zeros(2*N+2,nt);
m_true(N+2:2*N+2,1)=height_reg_weight; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT INITIAL WAVE HEIGHT %%%%%%%%%%%%%%%%%%%%%
fig=fig+1;
figure(fig)
%plot(x_reg,height_reg_weight,x_reg,height_reg);
plot(x_reg,height_reg_weight,'LineWidth',5);
%title('Initial wave height','Fontsize',20)
%legend('after weighting','before weighting') 
xlabel('Easting [km]','FontSize',15)
ylabel('Wave height [m]','FontSize',15) 

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUILD PROPAGATION MATRICES %%%%%%%%%%%%%%%%%%
T=compute_T_matrix(N,alpha,M,D,SAT,A,cdiss);
T=sparse(T);

%%%%%%%%%%%%%%%%%%%%%%%% PROPAGATION / WAVE MODELING %%%%%%%%%%%%%%%%%%%%%%
wantplot=1;
wait=0;
fig=fig+1;

for it=2:nt
    
    % update time
    t = t0 + (it-1) * dt;  
    
    % integrate one time step
    m_true(:,it)=T*m_true(:,it-1);    

    %%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%% 
    if (wantplot==1)
        if mod(it,50) == 0
            figure(fig)
            time_str = num2str(t);
            plot(x,m_true(N+2:2*N+2,it),'o')
            s = ['wave height at time: ',time_str,'(s)']; 
            title(s)
            ylim([-5 5])
            xlim([0 500])
            if (wait) 
               pause
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%% END PLOT %%%%%%%%%%%%%%%%     
end

fprintf('finished generating synthetic data \n')

%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA EXTRACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x_receivers_grid_temp,x_receivers_temp]=recv_init(nr,or,dr,dx,N); % set receivers' coordinates
H=compute_H_matrix(N,x_receivers_grid_temp); % compute data extraction matrix
H=sparse(H);
data=H*m_true; % extract data at receivers' locations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QC DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recv_coord_sim=[N+1+or:dr:N+1+or+dr*(nr-1)]';
% for it=1:nt
%     if mod(it,50) == 0
%         subplot(3,1,1)
%         plot(recv_coord,data(:,it))
%         subplot(3,1,2)
%         plot(recv_coord,m_true(N+1+or:N+1+or+dr*(nr-1),it))
%         subplot(3,1,3)
%         plot(recv_coord,m_true(N+1+or:N+1+or+dr*(nr-1),it)-data(:,it))
%         pause
%     end
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KALMAN FILTERING %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_variance=1.0e-2; % noise variance (square it for stand dev)
noise=sqrt(data_variance).*randn(nd,nt);
data=data+noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COVARIANCE MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%
% Q matrix (cov matrix of state evolution equation error)
Q_diag=1.0e-6;
Q=Q_diag*eye(nm,nm);

% R matrix (data covariance matrix)
R_diag=data_variance;
R=R_diag*eye(nd,nd);

% S_0 (initial model covariance matrix)
S_0_diag=0;
S_0=S_0_diag*eye(nm,nm);

% call Kalman operator function
fprintf('starting Kalman filter \n')
[m,S,G,PD,Nu]=kalman_filter_op(nt,nm,T,H,Q,R,S_0,data);
fprintf('finished Kalman filter \n')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('./simulations/s1_N','N')
save('./simulations/s1_nm1','nm1')
save('./simulations/s1_nm','nm')
save('./simulations/s1_nt','nt')
save('./simulations/s1_time','time')
save('./simulations/s1_x','x')
save('./simulations/s1_m','m')
save('./simulations/s1_m_true','m_true')
save('./simulations/s1_data','data')
save('./simulations/s1_x_receivers','x_receivers')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;

N=load('./simulations/s1_N');N=N.N;
nm1=load('./simulations/s1_nm1');nm1=nm1.nm1;
nm=load('./simulations/s1_nm');nm=nm.nm;
nt=load('./simulations/s1_nt');nt=nt.nt;
time=load('./simulations/s1_time');time=time.time;
x=load('./simulations/s1_x');x=x.x;
m=load('./simulations/s1_m');m=m.m;
m_true=load('./simulations/s1_m_true');m_true=m_true.m_true;
data=load('./simulations/s1_data');data=data.data;
x_receivers=load('./simulations/s1_x_receivers');x_receivers=x_receivers.x_receivers;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% WAVE HEIGHT AND FLUX %%%%%%%%%%%%%%%%%%%%%%%%%
wait=1;
for it=1:nt

    if (mod(it,10)==1)
           
        time(it)
        %figure(1)
        X1a=m(nm1+1:nm,1,it);
        X1b=m(nm1+1:nm,2,it);
        %X1residual=data(:,it)-m(nm1+1:nm,2,it);
        X1residual=m_true(nm1+1:nm,it)-m(nm1+1:nm,2,it);
        
%         %subplot(2,1,1)
%         figure(1)
          plot(x,X1b,x_receivers,data(:,it),x,m_true(nm1+1:nm,it));   
%         plot(x,X1b,x_receivers,data(:,it),x,m_true(nm1+1:nm,it),'LineWidth',4); 
%         %iter=num2str(it);
%         %s=['Wave height at iteration # ',iter];
%         %title(s,'Fontsize',20)
        legend('reconstructed height','data','true')    
        axis([0,400,-4,8])
        xlabel('Easting [km]','FontSize',15)
        ylabel('Wave height [m]','FontSize',15)     
        
        
        X2a=m(1:nm1,1,it);
        X2b=m(1:nm1,2,it);
        X2residual=m_true(1:nm1,it)-m(1:nm1,2,it);
    
        %subplot(2,1,2)
%         figure(2)
%         plot(x,X2b,x,m_true(1:nm1,it),x,X2residual);
%         s=['Flux at iteration # ',iter];
%         title(s,'Fontsize',20)  
%         legend('reconstructed flux','true flux','error')
%         axis([0,600,-1,1])
%         xlabel('Easting [km]','FontSize',20)
%         ylabel('Particle flux [10^3*m^2/s]','FontSize',20)     
        
%         time(it)
% 
%         figure(2)
%         plot(x,m_true(nm1+1:nm,it),'LineWidth',5);  
%         axis([0,600,-4,8])
%         xlabel('Easting [km]','FontSize',15)
%         ylabel('Wave height [m]','FontSize',15)    
%         
%         figure(3)
%         plot(x,m_true(1:nm1,it),'g','LineWidth',5);  
%         axis([0,600,-0.8,0.8])
%         xlabel('Easting [km]','FontSize',15)
%         ylabel('Particle flux [10^3*m^2/s]','FontSize',15)       
        
    
%         subplot(4,1,3)
%         plot(x,diag(S(1:nm1,1:nm1,2,it)));
%         s=['Height variance at iteration # ',iter];
%         title(s,'Fontsize',20)  
%         legend('reconstructed flux','true flux','error')
%         axis([0,600,-1,1])
%         xlabel('Easting [km]','FontSize',20)
%         ylabel('Particle flux [10^3*m^2/s]','FontSize',20)  
% 
%         subplot(4,1,4)
%         plot(x,diag(S(nm1+1:nm,nm1+1:nm,2,it)));
%         s=['Flux variance at iteration # ',iter];
%         title(s,'Fontsize',20)  
%         legend('reconstructed flux','true flux','error')
%         axis([0,600,-1e-3,1e-3])
%         xlabel('Easting [km]','FontSize',20)
%         ylabel('Particle flux [10^3*m^2/s]','FontSize',20)        
        
         if (wait)
             pause
         end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D PLOT OF TSUNAMI HEIGHT %%%%%%%%%%%%%%%%%%%%
ntplot=nt;
nxplot=nm1;
htrue=m_true(nm1+1:nm,1:ntplot)';
hrec=m(nm1+1:nm,2,1:ntplot);
hrec=reshape(hrec,nm1,ntplot)';
hdiff=htrue-hrec;
[x_mesh,t_mesh]=meshgrid(x(1:nxplot,1),time(1:ntplot));

% fig=fig+1;
% figure(fig)
%subplot(1,3,1)
figure(1)
xmax=300;
pcolor(x_mesh(:,1:xmax),t_mesh(:,1:xmax),htrue(:,1:xmax));
set(gca,'YDir','reverse')
%shading flat
shading interp
set(gca,'CLim',[-4,8]);
%h=colorbar;
%set(get(h,'title'),'string','Height [m]','Fontsize',15);
% title('True height','Fontsize',15)
%xlabel('Easting [km]','Fontsize',15)
%ylabel('Time [s]','Fontsize',15)

% fig=fig+1;
% figure(fig)
%subplot(1,3,2)
figure(2)
pcolor(x_mesh(:,1:xmax),t_mesh(:,1:xmax),hrec(:,1:xmax));
set(gca,'YDir','reverse')
%shading flat
shading interp
set(gca,'CLim',[-4,8]);
h=colorbar;
set(get(h,'title'),'string','Height [m]','Fontsize',15);
% title('Reconstructed height')
%xlabel('Easting [km]','Fontsize',15)
%ylabel('Time [s]','Fontsize',15)

% fig=fig+1;
% figure(fig)
% subplot(1,3,3)
figure(3)
pcolor(x_mesh(:,1:xmax),t_mesh(:,1:xmax),hdiff(:,1:xmax));
set(gca,'YDir','reverse')
%shading flat
shading interp
set(gca,'CLim',[-4,8]);
%h=colorbar;
%set(get(h,'title'),'string','Height [m]','Fontsize',15);
%title('Reconstruction error')
%xlabel('Easting [km]','Fontsize',15)
%ylabel('Time [s]','Fontsize',15)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D PLOT OF TSUNAMI FLUX %%%%%%%%%%%%%%%%%%%%%%
ntplot=nt;
nxplot=nm1;
ftrue=m_true(1:nm1,1:ntplot)';
frec=m(1:nm1,2,1:ntplot);
frec=reshape(frec,nm1,ntplot)';
fdiff=ftrue-frec;
[x_mesh,t_mesh]=meshgrid(x(1:nxplot,1),time(1:ntplot));

subplot(1,3,1)
pcolor(x_mesh,t_mesh,ftrue);
set(gca,'YDir','reverse')
shading flat
shading interp
set(gca,'CLim',[-1,1]);
h=colorbar;
set(get(h,'title'),'string','Wave Height [m]','Fontsize',20);
title('True flux')
xlabel('Easting [km]','Fontsize',20)
ylabel('Time [s]','Fontsize',20)

%fig=fig+1;
%figure(fig)
subplot(1,3,2)
pcolor(x_mesh,t_mesh,frec)
set(gca,'YDir','reverse')
shading flat
shading interp
set(gca,'CLim',[-1,1]);
h=colorbar;
set(get(h,'title'),'string','Wave Height [m]','Fontsize',20);
title('Reconstructed flux')
xlabel('Easting [km]','Fontsize',20)
ylabel('Time [s]','Fontsize',20)

%fig=fig+1;
%figure(fig)
subplot(1,3,3)
pcolor(x_mesh,t_mesh,fdiff)
set(gca,'YDir','reverse')
shading flat
shading interp
set(gca,'CLim',[-1,1]);
h=colorbar;
set(get(h,'title'),'string','Wave Height [m]','Fontsize',20);
title('Flux reconstruction error')
xlabel('Easting [km]','Fontsize',20)
ylabel('Time [s]','Fontsize',20)