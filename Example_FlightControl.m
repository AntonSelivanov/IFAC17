% This MATLAB program checks the feasibility of LMIs from Theorem 1 and Remark 5 of the paper 
% A. Selivanov, E. Fridman, and A. Fradkov, "Event-triggered adaptive control of minimum-phase systems," in 20th IFAC World Congress, 2017. 

%% Parameters of the plane 
a1=[-1.5 -.7]; % uncertain
a2=33; 
a3=-1.3; 
b1=19/15; 
b2=19; 
kb=-1.5e-4; 
Tb=1/65; 
xib=.01; 

%% Plant
for i=1:length(a1) 
    Ap{i}=[a1(i) 1 0 0 0; 
           a2 a3 0 0 0; 
           0 1 0 0 0; 
           0 0 0 -2*xib/Tb 1;
           0 0 0 -1/Tb^2 0]; %#ok<*SAGROW>
end
Bp=[b1; b2; 0; 0; kb/Tb^2]; 
Cp=[0 0 1 1 0]; 
gp=1; 
r=2; % relative degree

%% Shunt 
As=-14; 
Bs=2; 
Cs=1; 

%% Augmented system
for i=1:length(Ap)
    A{i}=blkdiag(Ap{i},As); 
end
B=[Bp; Bs]; 
C=blkdiag(Cp,Cs); 
g=[gp; 1]; 

%% Tuning parameters 
sigma=4e-3; % Event-triggering threshold
M=35; 
kStar=30; 

%% Periodic sampling (Theorem 1 with sigma=0)
h0=.005; 
hmax=fminsearch(@(h) -LMI_IFAC16_th1(A,B,C,g,r,kStar,h,0,M),h0);
if hmax==h0
    display('Periodic sampling: not feasible'); 
else
    display(['Periodic sampling: h=' num2str(hmax)]); 
end

%% Switching event-triggering (Theorem 1)
h0=.005; 
hmax=fminsearch(@(h) -LMI_IFAC16_th1(A,B,C,g,r,kStar,h,sigma,M),h0);
if hmax==h0
    display('Swithcing event-triggering: not feasible'); 
else
    display(['Switching event-triggering: h=' num2str(hmax)]); 
    [~,P]=LMI_IFAC16_th1(A,B,C,g,r,kStar,hmax,sigma,M); 
end

%% Periodic event-triggering (Remark 5)
h0=.005; 
hmax=fminsearch(@(h) -LMI_IFAC16_rem5(A,B,C,g,r,kStar,h,sigma,M),h0);
if hmax==h0
    display('Periodic event-triggering: not feasible'); 
else
    display(['Periodic event-triggering: h=' num2str(hmax)]); 
    [~,P]=LMI_IFAC16_rem5(A,B,C,g,r,kStar,hmax,sigma,M); 
end
