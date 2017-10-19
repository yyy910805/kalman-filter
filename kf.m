function [ wa, Va ] = kf ( All_Verbose, Tfinal, CN, do_filter )
%KF   Kalman filter code for linear advection equation
%
%       Usage:  kf ( All_Verbose, Tfinal, CN, do_filer )
%
%       Inputs:
%
%       All_Verbose   -  when diff than zero, will prompt a full set of
%                        large set of options
%       Tfinal        -  total intergration time
%       CN            -  Courant number
%       do_filter     -  0=carry estimation using KF; 
%                        1=no estimation
%
%       List of questions user will be prompted to: 
%
%                                            Default values
%                                            --------------
%       Total time of integration        -     1 time unit
%
%       Courant number                   -       1
%       Initial error standard deviation -       0.05
%
%       Parameters related to Obs Pattern:
%       Observation frequency            -       dt  (every model time step)
%       Observation error standard dev.  -       0.02
%       Observation sparsity             -       every grid point
%
%       Plotting frequency               -       0 (only in the beginning
%                                                   and end of a run)
%                                                -  set this to a multiple
%                                                   of dt and you get plots
%                                                   with this frequency
%
%  25feb98  R. Todling   Initial code.

figure(1), clf
if ( do_filter == 0 ), figure(2), clf, end
echo off
%
%  Define length of assimilation
%
% Tfinal = input('Total time of assimilation experiment (1 time unit):');
% if isempty(Tfinal), Tfinal=1, end
%
%  Define parameters
%
dt   = 0.05;        dx   = 0.20;% 0.05 %0.20;
La   = -2;          Lb   = 2;
x    = La:dx:Lb-dx;

jdim = size(x,2);
wt   = zeros(jdim,1);
%
%  Define advection-diffusion eq parameters
%
% CN = input('Enter Courant number (default=1):');
% if isempty(CN), CN=1, end
echo on
U  = CN * dx / dt; %  Advetion speed
MU = 0;            %  Diffusion coefficient
echo off
global Courant_Number Diffusivity
Courant_Number  = CN;  % U  * dt / dx
Diffusivity     = MU * dt / (dx*dx);
%Re = Courant_Number / Diffusivity

%
% Build initial condition
%
w0 = sqrwv ( x, La, Lb );
%
%  Initialize error covariances
%
ID = eye(jdim); ZERO = zeros(jdim); 
Ld = 0.5*abs(La);
%
%  Define a correlation matrix and its Cholesky decomposition
%
Cor0 = gcorr ( 'gauss', Lb-La, Ld, jdim, 0 );
[V,D]=eig(Cor0);  sCor0 = V * sqrt(D); 
%
%  Set the true state to the analytic function w0 and define
%  initial error magnitude
%
wt = w0;
stdAd = 0.05; stdA=stdAd;
if (All_Verbose~=0)
    stdA= input('Enter initial error magnitude (0.05):');
end

if isempty(stdA)
    stdA=stdAd;
end
pa = stdA^2 * Cor0;
wa = zeros(jdim,1);

if (All_Verbose~=0), Iwhich = input('Case to study: 1=wave generation; 2= general:');
   if isempty(Iwhich), Iwhich=1; end
   if (Iwhich == 1), wa = zeros(jdim,1); end
   if (Iwhich == 2), wa = w0 + stdA * sCor0 * randn(jdim,1); end
end
%
%  Because of linearity, we can save the dynamics in a matrix
%
psi = getpsi('upwind',jdim);
%
%  Set specifics about model error
%
stdQd = 0; stdQ=stdQd;
if (All_Verbose~=0), stdQ= input('Enter model error magnitude(std dev, set to null):');
   if isempty(stdQ), stdQ=stdQd; end
end
qq = stdQ^2 * Cor0;

%
%  Define observing network info
%
tobs = 99999;
if ( do_filter == 0 )
   tobs  = input('Enter observation frequency (this*dt):');
   
   if isempty(tobs), tobs=1; end 
   tobs   = tobs * dt;
   nobs  = input('Enter observation sparsity (obs over left-half):');
   
   if isempty(nobs), nobs=-jdim/2; end
   stdO = input('Enter observation error std dev (0.02):');
   
   if isempty(stdO), stdO = 0.02; end
%
%  Construct observation error covariance and observation matrix
%
	[hh,io] = obspat(nobs,jdim);
	nobs    = length(io);
	rr      = stdO^2 * eye(nobs);
   
end   
%
%  Gather I/O information
%
Plot_Nowd = 999999; Plot_Now=Plot_Nowd;
if (All_Verbose~=0)
    Plot_Now  = input('Enter plotting frequency (this*dt)');
end
if isempty(Plot_Now)
    Plot_Now = Plot_Nowd;
end
 
if (Plot_Now ~= 999999), figure(3), clf, end
%
%  Now plot some initial stuff: initial truth, initial estimate,
%      observation pattern
%        
ymaxT = max(wt); ymaxA = max(wa); ymax=max(ymaxT,ymaxA); 
yminT = min(wt); yminA = min(wa); ymin=min(yminT,yminA);
figure(1),subplot(211), hold on,plot(x,wt,'c-'),axis([x(1) x(jdim) ymin ymax]) ;
figure(1),subplot(211), hold on,plot(x,wa,'b-');
if(do_filter==0)
    figure(1),subplot(211), hold on,plot(x(io),zeros(length(io),1),'y*');
    figure(1),subplot(211), hold on, ...
    legend('Truth at T0','Estimate at T0','Obs Pattern'), hold off;	  
else
    legend('Truth at T0','Estimate at T0'), hold off;	  
end       
%echo on;
%pause % Strike any key for plot.
%echo off;
%  
%  Start Main Time Loop
%
kcnt = 1;
Va   = sqrt(sum(diag(pa)));
for t = dt:dt:Tfinal
    
   if ( mod(t,Plot_Now*dt) == 0 )
      figure(3),plot(x,wt,'r-'), hold on;
      figure(3),plot(x,wa,'g-'), hold on;
      %echo on;
      pause % Strike any key for plot.
      %echo off;
   end
% 
%  Evolve true state
%
    wt = upwind ( wt, t )  +  stdQ * sCor0 * randn(jdim,1);
%
%  Evolve estimate 
%
    wf = upwind ( wa, t );
%
%  Evolve covariance
%
    pf  = psi * pa * psi'  +  qq;
    %figure(3), plot(x,pf(:,1)), hold on, pause;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Observation/Assimilation time
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ( mod(t,tobs) == 0)
      kcnt = kcnt + 1; 
%
%  Generate observations by adding random error
%
      wo = hh * wt  +  stdO * randn(nobs,1);
%
%  Compute Kalman Gain
%
      kk = pf * hh' * inv( hh*pf*hh' + rr );
%
%  Update State and Error Covariance 
%
      wa =  wf  +  kk * ( wo - hh*wf );
      pa = ( ID - kk*hh ) * pf * ( ID - kk*hh )'  +  kk * rr * kk';
      
  else

      wa = wf;
      pa = pf;

  end  %  End Observation loop
  Va = [ Va, sqrt(sum(diag(pa)))];

end   %  End Main Time Loop
%
%  Always plot the file results
%
ymaxT = max(wt); ymaxA = max(wa); ymaxN=max(ymaxT,ymaxA); ymax=max(ymax,ymaxN);
yminT = min(wt); yminA = min(wa); yminN=min(yminT,yminA); ymin=min(ymin,yminN);
figure(1),axis([x(1) x(jdim) ymin ymax])
figure(1),subplot(212), hold on,plot(x,wt,'r-'), ...
          axis([x(1) x(jdim) ymin ymax])
figure(1),subplot(212), hold on,plot(x,wa,'g-')
figure(1),subplot(212), hold on, ...
   legend('Truth at Tf','Estimate at Tf'), hold off;
if (do_filter==0)
    figure(2),plot(Va);title('Expected RMS');
end
