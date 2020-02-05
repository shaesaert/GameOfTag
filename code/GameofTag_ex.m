%% \eps,\delta ROBUST CORRECT_BY_DESIGN CONTROL 
% IN THIS WORK, WE GIVE SOME EXPERIMENTAL ROUTINES FOR THE COMPUTATION OF
% \EPS \DELTA ROBUST SPECIFICATIONS.
% Copyleft @ Sofie Haesaert 10/2017 updated 2018-2019
clear all
close all
clc

%% DEFINE THE SPECIFICATION OF INTEREST
% we currently focus on the following specification 
% F G^n2 K
% K = is a polytope 
n2=3;
DFA.S= 1:n2+1; % set of states
DFA.S0 = [1]; % initial state
DFA.Act = ['k','nk'];
DFA.F = DFA.S(end); % target states
DFA.Trans=[2:DFA.F,DFA.F;ones(1,length(2:DFA.F)),DFA.F]'; % column = control , row = state, value = next state
% Can be verified with 
% https://spot.lrde.epita.fr/trans.html
% 

%% DEFINE THE STOCHASTIC MODEL
% The original/concrete model is an LTI model with dimensions m,mw,n and parameters
% a,b,c and matrices A,B,Bw and C. 
m=1;
mw=1;
a=.3;
b=.8;
c=.8;
n=3;
A = [1 -a a
    0 b  0
    0 0   c];
eig(A)
B = [-a*.1;1;0];
Bw=[a*.02;0;.1];
C =[1 0 0];


LTI_concrete = ss(A,[B Bw],C,[],1); % make it a dt state space
K_pol = Polyhedron([-2,2]'); % Ideal following distance (target) 
             

%% GO FROM 3 DIMENSIONAL TO 1 DIMENSIONAL LTI MODEL
% Reduce 3 dimensional model to 1 dimensional LTI model
opt = balredOptions('Offset',.001');  
[Ml,~, K]=dare(LTI_concrete.A,LTI_concrete.B(:,1),.5*eye(n)... This term enforces a whicht on the precision
    ,0.08)  ;
LTI_concrete_cl = ss(A-B*K,[B Bw],C,[],1); % make it a dt state space

[sysred] = balred(LTI_concrete_cl,1); 

if sysred.c(end)~= 1
    disp(' sysred.c(end)~= 1')
    disp('  sysred.c')
    sysred.c;
    
    T=eye(length(sysred.c));
    T(length(sysred.c),length(sysred.c))= sysred.c(end);
    sysred=ss2ss(sysred,T);
end
    sysred.d=zeros(1,m+mw);
     
 
LTI_abstract.A= sysred.A;
LTI_abstract.B= sysred.B(:,1);
LTI_abstract.Bw= sysred.B(:,2:end); 
LTI_abstract.C= sysred.C;


%% DETERMINE THE characteristics for the gridding
LTI_abstract.U = Polyhedron([-.3,.3]);
LTI_abstract.X = Polyhedron([-10,10]);  
% - X = Polytope of the to be gridded part of the state space
% - U = Polytope of the to be gridded part of the input space

rad=.05; % The size of the grids
% prelim evaluation of the precision
[eps,del,~,~,~,~ ] = epsdel_compute(LTI_concrete,LTI_abstract, Ml,K,rad,.01);

fprintf('Initial epsilon =  %d\n delta = %d, \n',eps,del);


%%  GO FROM 1 DIMENSIONAL LTI MODEL TO MDP
% BASED ON THE GIVEN GRIDDING PRECISION GRID THE LTI MODEL
% for this we call the gridding function

nu=20; % number of grid points for the input space

[MDP,rad] = gridding(LTI_abstract,  rad, nu);
% with 
% - LTI_abstract = the abstract 1 dimensionsal model
% - diam = the gridding diameter (euclidean norm for higher dimensions)
% - nu = the number of grid points/representative points in the input space

[eps,del,Q,R,P,M ]= epsdel_compute(LTI_concrete,LTI_abstract, Ml,K,rad,.03);

fprintf('Final epsilon =  %d\n delta = %d, \n',eps,del);


% Compute Phat to prject initial states
Phat= (P'*M*P)\P'*M;


%% GO FROM DFA TO NFA, WHOSE TRANSITIONS TAKE AS INPUT THE STATE OF THE MDP
% Now, we need a DFA that represents this. 
% the NFA contains states, initial states, goal states,  transitions
% relations
NFA=  NFA_eps(DFA,eps,MDP,K_pol);




%% FOR THE 1 DIMENSIONAL MODEL COMPUTE THE \DELTA-ROBUST SATISFACTION AND THE POLICY
%  
[p,mu] = del_reach(MDP, NFA, del);

plot(MDP.z_rep,p)
title('Robust satisfaction probability')




%% Do a simulation with the refined controller
figure
for run=1:10 % 10 runs
%1. Initiate
x_2=[2.45;2.5;1.3]; % initial state of orig/concrete model 
x_1=Phat*x_2; % reduced order model initial state
q=DFA.S0; % initiate the DFA
N=10; % simulation horizon 
% start simulation:
for t =1:N
     % update compute q based on y_2
    q(t+1)=DFA.Trans(q(t),1)*K_pol.contains(C*x_2(:,t)) +DFA.Trans(q(t),1)*(~K_pol.contains(C*x_2(:,t))  );
    
    % find representative state in abstract model 
    [v,maxrep] =min(abs(MDP.z_rep-x_1(:,t)*ones(1,length(MDP.z_rep))));
    % abstract input
    u1=mu(q(t+1),maxrep);
    % map back to continuous reduced order state 
    x_1c(:,t) =MDP.z_rep(maxrep);
    
    % concrete input
    u2 = R*u1+Q*x_1c(:,t)-K*(x_2(:,t)-P*x_1(:,t));
    
    % noise realization that drives state evolution
    w=randn(1,1);
    x_1(:,t+1)=LTI_abstract.A*x_1(:,t)+LTI_abstract.B*u1...
                +LTI_abstract.Bw*w;
    
    x_2(:,t+1)=A*x_2(:,t)+B*u2+Bw*w;
    % Remark: 
    % for control synthesis w has to be computed from x_1(t) and x_1(t+1), that is, 
    % w = (LTI_abstract.Bw'*LTI_abstract.Bw)^(-1)*LTI_abstract.Bw'*... 
    % (x_1(:,t+1) - LTI_abstract.A*x_1(:,t)+LTI_abstract.B*u1)
    
    % map x_1 to representative grid point
    [~,maxrep] =min(MDP.z_rep-x_1(:,t+1)*ones(length(MDP.z_rep)));
    
end

% Plot simulation:

plot(1:N,x_1c(1,:),'x')
hold on
plot(1:N+1,x_2(1,:))

end