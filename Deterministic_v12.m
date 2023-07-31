%% Main_Alg
close all;
clear all;
clc;
format short
diary off
%
for run_iter= 2
    ins_iter= 5;
    MatName= ['ins_', num2str(run_iter), '_iter_', num2str(ins_iter), '_V8'];
    if run_iter==4 || run_iter==7
        detergap= .1;
    elseif run_iter==10
        detergap= .15;
    else
        detergap= 10^(-4);
    end
    DiaryName= [MatName, '_Deter_v12_', num2str(detergap),'.txt'];
%     diary(DiaryName)
    disp(DiaryName);
    disp(['Cplex Ver: ',num2str(getVersion(Cplex))])
    CLOCK=clock;
    disp(['Clock: ',num2str(CLOCK(4)),':',num2str(CLOCK(5))]);
    disp('--------------------------------');
    %% Main Inputs
    %% Algorithm tuning
    load([pwd, '\inputs\ins_', num2str(run_iter), '\', MatName]);  % retailer sort: {1,1 , 2,2 , 3,3}, time sort: {1,2 , 1,2 , 1,2}
    N=mat.N; T=mat.T; D=T*N; K=mat.K; M=mat.M;
    timelimit=400*D;
    %% Sets
    V=M+N;
    disp(['N: ',num2str(N)])
    disp(['T: ',num2str(T)])
    disp(['K: ',num2str(K)])
    disp(['M: ',num2str(M)])
    %% Variavles
    % Initial inventory level of each brewery shuold be lower than total filled bottles demand
    purch= sdpvar(M,1);
    m= sdpvar(M,T,'full');
    r= sdpvar(M,T,'full');
    Iful= sdpvar(M,T+1,'full');        %Iful(:,1) related to time zero
    Iemp= sdpvar(V,T+1,'full');        %%Iemp(:,1) related to time zero
    f= sdpvar(N,K,T,'full');
    e= sdpvar(N,K,T,'full');
    v= sdpvar(V,V,K,T,'full');
    u= sdpvar(V,V,K,T,'full');
    x= binvar(V,V,K,T,'full');
    %% Parameters
    delta= [8;0];%mat.deter_delta;
    delta=[zeros(M,T);delta];
    pi= [8;0];%mat.deter_pi;
    pi=[zeros(M,T);pi];
    Q= mat.Q;              % Capacity of vehicle
    FilCap= mat.FilCap;    % Filling capacity
    Cpur= mat.Cpur;        % purchasing cost
    Chyg= mat.Chyg;        % hygienization cost
    Cins= mat.Cins;        % inspection cost
    Cdis= mat.Cdis;        % disposal reimburse
    FixC=mat.FixC;         % fix cost of vehicles
    C= mat.C;
    He= mat.He;
    Hf= mat.Hf;            % Inventory holding cost of filled bottles
    Lf= mat.Lf;            % Storage capacity for filled bottles
    Le= mat.Le;            % Storage capacity for empty bottles
    alpha= mat.alpha;      % Hygienization rate
    Iful(1:M,1)= mat.Iful;
    Iemp(1:V,1)= mat.Iemp;
    %% constraints
    Constraints = [
        sum(sum(x(1:M,1:M,:,:),3),4) == 0 ;    % Shipment among Breweries is not allowed
        diag(sum(sum(x(M+1:V,M+1:V,:,:),3),4)) == 0 ;    % Routing from a retailer to itself is not allowed
        r <= alpha .* Iemp(1:M,1:T);%Cons 6
        m + r <= repmat(FilCap,1,T); %Cons 7
        Iful(1:M,2:T+1) <= repmat(Lf,1,T);%Cons 9
        Iemp(1:V,2:T+1) <= repmat(Le,1,T);%Cons 10
        sum(x(1:M,M+1:V,1:K,1:T),[1,2]) <= 1;%Cons 11
        sum(sum(x(M+1:V,:,:,1:T),2),3) <= 1;%Cons 12
        sum(x,2) - permute(sum(x,1),[2 1 3 4]) == 0;%Cons 13
        v + u <= Q .* x;%Cons 16
        purch >= 0; m >= 0; r >= 0; Iful(:,2:T+1) >=0;%Cons 17
        Iemp(1:V,2:T+1) >=0;%Cons 18
        f >= 0; e >= 0;%Cons 19
        v >= 0; u >= 0 ];%Cons 20
    for period=1:T
        Constraints = [Constraints;
            Iful(1:M,period) + m(1:M,period) + r(1:M,period) - sum(sum(v(1:M,M+1:V,:,period),2),3) == Iful(1:M,period+1);%Cons 2
            Iemp(1:M,period) + permute(sum(sum(u(M+1:V,1:M,:,period),1),3),[2 1 3 4]) - (r(1:M,period)/alpha)   == Iemp(1:M,period+1); %Cons 3
            Iemp(M+1:V,period) - sum(e(:,:,period),2) + pi(M+1:V ,period)   <= Iemp(M+1:V,period+1);%Cons 4
            sum(f(:,:,period),2) >= delta(M+1:V ,period);%Cons 5
            m(:,period) <= purch- sum(m(:,1:period-1),2)]; %Cons 8
    end
    for k=1:K
        Constraints = [Constraints;
            permute(sum(v(:,M+1:V,k,1:T),1),[2 1 3 4]) - sum(v(M+1:V,:,k,1:T),2) == f(1:N,k,1:T);%Cons 14
            sum(u(M+1:V,:,k,1:T),2) -  permute(sum(u(:,M+1:V,k,1:T),1),[2 1 3 4]) == e(1:N,k,1:T)];%Cons 15
    end
    %% add extra
    for period=1:T
        Constraints = [Constraints;
            sum(f(:, :, period), 1)' <= permute(sum( Q(1:M, 1:V, :, period) .* x(1:M, 1:V, :, period), [1 2]),[3,1,2]);   %28_paper
            sum(e(:, :, period), 1)' <= permute(sum( Q(1:M, 1:V, :, period) .* x(1:M, 1:V, :, period), [1 2]),[3,1,2]) ];  %30_paper
        for i=M+1:V
            Constraints = [Constraints;
                sum(x(i, 1:V, :, period),2) <= sum( x(1:M, 1:V, :, period), [1 2]) ];  %31_paper
        end
    end
    maxQ= max(Q(1,1,:,1));
    for period=1:T
        Constraints = [Constraints;
%             sum(sum(sum(x(M+1:V, 1:V, :, 1:period),2), 3), 4) >= ceil( sum(delta(M+1:V,1:period),2) ./ maxQ );  %25_paper (remove: because in our study, a retailer can be visited just by one vehicle, but in the base study (Qiu 2018) there is no limitation)
            sum(x(1:M, 1:V, :, 1:period),'all') >= ceil(sum(delta(M+1:V,1:period),[1,2]) ./ maxQ ) ];  %26_paper
%             sum(sum(sum(x(M+1:V, 1:V, :, 1:period),2), 3), 4) >= ceil(( Iemp(M+1:V,1) + sum(pi(M+1:V,1:period),2) - Le(M+1:V) ) ./ min(maxQ,Le(M+1:V)) ) ];  %27_paper (remove: because in our study, a retailer can be visited just by one vehicle, but in the base study (Qiu 2018) there is no limitation)
        for s=0:period-1
            Constraints = [Constraints;
                sum(delta(M+1:V, period-s:period),2) .* (1- sum(sum(sum(x(M+1:V, :, :, 1:period), 2), 3), 4)) <= 0;  %33_paper
                Le(M+1:V) - sum(pi(M+1:V, period-s:period),2) .* (1- sum(sum(sum(x(M+1:V, :, :, 1:period), 2), 3), 4)) >= Iemp(M+1:V,period-s) ];  %34_paper
        end
    end
    %% ObjFunc
    Objective= sum(sum(Cpur * purch)) + sum(sum(Chyg * r)) ...
        + sum(sum((Cins-((1-alpha) * Cdis)) * Iemp(1:M,1:T))) ...
        + sum(sum(sum(sum(C .* (v + u))))) ...
        + sum(sum(sum(sum(FixC .* x)))) ...
        + (sum(He' * Iemp(:,2:T+1))) ...
        + (sum(Hf' * Iful(:,2:T+1))) ;
    for tt=1:T
        Objective=Objective+ He(1:M)'*(purch- sum(m(:,1:tt),2));
    end
    timing=tic;
    sol=optimize(Constraints, Objective, sdpsettings('solver','cplex','verbose',3 ,'cplex.timelimit', timelimit,'cplex.mip.tolerances.mipgap', detergap ));
    disp(' ')
    disp(' ')
    disp(['Deterministic time= ', num2str(sol.solvertime)])
    disp(['sol_status= ', num2str(sol.problem)])
    total_time=toc(timing);
    %% Result and Figures
    Obj_sol= value(Objective);
    if sol.problem ~= 0
        Obj_sol= Obj_sol + Inf;
    end
%     TIME=clock;
%     disp(['Clock: ',num2str(TIME(4)),':',num2str(TIME(5))]);
    %
    val.purch= value(purch);
    val.m= value(m);
    val.r= value(r);
    val.Iful= value(Iful);
    val.Iemp= value(Iemp);
    val.f= value(f);
    val.e= value(e);
    val.v= value(v);
    val.u= value(u);
    val.x= value(x);
    val.delta=delta;
    val.pi=pi;
    %
    disp(['Solution value= ', num2str(Obj_sol)]);
%     save([MatName,'_Deter_v12_',num2str(detergap),'.mat'],'val');    
%     zip([MatName,'_Deter_v12_',num2str(detergap)],{[MatName,'_Deter_v12_',num2str(detergap),'.mat'],[MatName, '_Deter_v12_', num2str(detergap),'.txt']});    
    diary off
end