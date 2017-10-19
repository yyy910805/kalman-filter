%SOLVE02  Solution of problem 2 in ...

%
%   Make selected plots of pde solution
%
clf reset
clc
echo on
%  This solves exercise 2 step by step. 
%  This is not completely self explanatory, you are expected to at least
%  have read the statement of the problem being solved. Basically, this
%  program makes a series of ordered calls to the program "kf.m", which
%  is provided together with this driver, to guide you through the solution
%  of each item in the problem.

pause % Hit any key to continue
clc
%
%  Simulation experiments:  In this part you are asked to compute the 
%                           numerical solution of the advection equation
%                           on a specified domain, with given parameters,
%                           for three different Courant numbers.

pause % Hit any key to see result for Courant number = 1
      % This corresponds to the EXACT solution after one revolution

echo off
clc

Tfinal    = 1;   % Final time of simulation
do_adv    = 1;   % Case of state evolution only (no estimation)
echo on

kf(0,Tfinal, 1, do_adv);   % Invoke simulation program 

% The top plot corresponds to the initial condition.
% The bottom plot corresponds to the state at the final time.
% The combination of the time step interval and the Courant
% number are such that the wave travels exactly one revolution,
% thus coming back to its initial position.

pause % Hit any key to see Result for Courant number = 0.95

kf(0,Tfinal, 0.95, do_adv);

% The solution differs from the previous one because the Courant number
% is different. Since the time step is kept as in the previsous case, and
% the speed is smaller than before, the wave now does not get to complete 
% one full revolution. Courant numbers smaller than one also result in 
% artificial (numerical) dissipation, causing to solution to spread.

pause % Hit any key to see Result for Courant number = 0.25

kf(0,Tfinal, 0.25, do_adv);

% In this case the advection speed is also different and 
% the artificial dissipation is quite strong.

pause % Hit any key to continue
clc

% The examples above are meant to be a quick demonstration of how much the 
% solution of the advection equation depends on the Courant number, and the 
% effect of that in generating artificial (numerical) dissipation.
 
%             ..... %  ......  %  ....... %  ........

% Before beginning the estimation experiments we need to construct a meaningful
% covariance matrix for the problem in question. In this case, we use a 
% homogeneous and isotropic formulation, representing Gaussianly distributed 
% errors.
%
% To do that in the domain of interest, it is important to make sure the matrix
% created numerically from the analytic form is indeed a covariance matrix. 
% It suffices to examine the spectrum of this covariance matrix, to make sure 
% all of the matrix eigenvalues are positive, thus assuring the matrix to be
% positive definite.

pause  % Hit any key to begin the solving the estimation problem
clc

x=-2:0.2:2-0.2; jdim = length(x);   % Define the domain (as for the PDE used before)

% Build a Gaussian covariance matrix

% The spectrum of Q shows that all eigenvalues are indeed positive, thus Q is
% a covariance matrix. The duplicity of eigenvalues is due to the periodicity 
% of the domain of interest.

pause  % Hit any key to continue
close all
clc
% 
%   Solution of the estimation problem using the Kalman filter
%   ==========================================================
%
%   The first case of estimating the state of the system when observing
%   over all grid points can be solved simply by hitting the return key
%   to three questions you will be prompted to below. This case correspond
%   to the default variables in the program "kf.m".
%   The three questions will be about:
%       (i) the observation frequency   - set to every time step
%      (ii) the observation sparsity    - set to one (obs over left-half)
%     (iii) the observation error level - set to a reasonable values
%

%  The half-observation coverage case:   In this case the initial estimate
%                                        is set to zero everwhere and the 
%                                        observation pattern is composed of
%                                        the grid points over the left-half 
%                                        the domain.
%

pause % Hit any key to continue
clc
                           % YOU ARE EXPECTED TO HIT THE RETURN KEY 
kf(0, Tfinal, 0.95, 0 );   % FOR ALL THREE FOLLOWING QUESTIONS

%  Since this is a linear problem, wasn't for the estimation (assimilation)
%  procedure the estimate would be zero at all times, since it had been set
%  to zero at the initial time.
%  It is the presence of the observations, combined with the assimilation 
%  machinery that allows for the estimate to become non-zero, and very close
%  to the true state.
%  The top panel in Fig. 1, shows the initial true state as well as the initial
%  guess (estimate) - set to zero. The bottom panel of Fig. 1 shows these
%  quantities at the final time of the assimilation period, where by then the
%  estimated solution is indeed quite close to the true state. Experiments
%  like this, in which the initial estimate is set to zero, are sometimes
%  referred to as "wave generation experiments", for obvious reasons.
%  Figure 2 shows the time evolution of the error standard deviation, extracted
%  from the diagonal of the filter error covariance matrix. Initially, this error
%  is quite high, but it becomes quite small at the last step, which is consistent
%  with the bottom plot in Fig. 1.

pause  % Hit any key to continue
clc

                           % In the next case, you want to increase the 
			   % observation level to 0.1:
			   %     HIT ENTER for the next 2 prompts and
			   %     ENTER 0.1 at the third prompt 
kf(0, Tfinal, 0.95, 0 );   

pause  % Hit any key to continue
clc

                           % To have both the observation level and its 
			   % frequency changed according to the statement 
			   % of the problem, ENTER 5 at the first prompt
			   %                 HIT RETURN at the second prompt
			   %                 ENTER 0.1 at the third prompt
kf(0, Tfinal, 0.95, 0 );
