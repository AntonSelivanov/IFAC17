function [hfeas,Pval]=LMI_IFAC16_rem5(A,B,C,g,r,kStar,h,sigma,M)
% This MATLAB program checks the feasibility of LMIs from Remark 5 of the paper 
% A. Selivanov, E. Fridman, and A. Fradkov, "Event-triggered adaptive control of minimum-phase systems," in 20th IFAC World Congress, 2017. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% A             - is a cell array of vetrices A_i=diag(A_{p,i},A_s)
% B, C          - parameters of the augmented system (6) 
% g             = col{g_p,1}
% r             - relative degree
% kStar         - large enough scalar such that (2) is feasible 
% h, sigma      - the event-triggering parameters from (11) 
% M             - positive scalar such that |k(0)-kStar|<M 

% Output: 
% hfeas         - equals h if feasible, 0 otherwise. Such output is convenient for fminsearch(); 
% Pval          - value of the decision matrix P

%% Notations 
if ~iscell(A)
    A={A}; 
end
n=size(A{1},1)-r+1; 
Cp=C(1:end-1,1:n); 
gp=g(1:end-1); 

%% Decision variables 
P=sdpvar(n+r-1); 
mu=sdpvar; 
w=sdpvar; 

%% LMIs 
LMIs=[P*B==C'*g,P>=0,mu>=0,w>=0]; 
for i=1:length(A) % loop over polytope vertices
    for a=[-M,M]
        Xi=blkvar; 
        Xi(1,1)=P*(A{i}-B*kStar*g'*C)+(A{i}-B*kStar*g'*C)'*P+mu*sigma*blkdiag((gp'*Cp)'*(gp'*Cp),zeros(r-1)); 
        Xi(1,2)=-P*B*kStar+C'*g*a; 
        Xi(1,3)=-P*B*kStar+C'*g*a; 
        Xi(1,4)=h*w*(A{i}'-C'*g*B'*(kStar+a))*[Cp'*gp;zeros(r-1,1)]; 
        Xi(2,2)=2*a-pi^2*w/4; 
        Xi(2,3)=2*a; 
        Xi(2,4)=-h*w*B'*(kStar+a)*[Cp'*gp;zeros(r-1,1)]; 
        Xi(3,3)=2*a-mu; 
        Xi(3,4)=-h*w*B'*(kStar+a)*[Cp'*gp;zeros(r-1,1)]; 
        Xi(4,4)=-w; 
        Xi=sdpvar(Xi); 
        
        LMIs=[LMIs, Xi<=0];  %#ok<*AGROW> 
    end    
end

%% Solution of LMIs
options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,[],options); 

hfeas=0; Pval=[]; 
if sol.problem == 0
    primal=check(LMIs); 
    if min(primal(2:end))>0
        Pval=value(P); 
        hfeas=h; 
    end
else
    yalmiperror(sol.problem) 
end
