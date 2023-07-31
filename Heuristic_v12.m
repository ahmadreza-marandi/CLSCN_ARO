%% Main_Alg
close all;
clear all;
clc;
format short
diary off
%
for run_iter= [3:7,9:10]
    for ins_iter= 1:5
        MatName= ['ins_', num2str(run_iter), '_iter_', num2str(ins_iter), '_V8'];
        if run_iter==4 || run_iter==7
            detergap= .1;
        elseif run_iter==10
            detergap= .15;
        else
            detergap= 10^(-4);
        end
        DiaryName= [MatName, '_H_v12_', num2str(detergap),'.txt'];
        diary([pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Heuristic_v12\', DiaryName])
        disp(DiaryName);
        disp(['Cplex Ver: ',num2str(getVersion(Cplex))])
        CLOCK=clock;
        disp(['Clock: ',num2str(CLOCK(4)),':',num2str(CLOCK(5))]);
        disp('--------------------------------');
        %% Main Inputs
        %% Algorithm tuning
        load([pwd, '\inputs\ins_', num2str(run_iter), '\', MatName]);
        N=mat.N; T=mat.T; D=T*N; K=mat.K; M=mat.M;
        %
        if D<=6
            dim= D;
        elseif D>6 && D<10
            dim= ceil(.5*D);
        else
            dim= 3;
        end
        if D<=4
            dim_scenario= D;
        elseif D>4 && D<10
            dim_scenario= ceil(.5*D);
        else
            dim_scenario= 3;
        end
        termint_time=400*D;    % On solution time
        siz_slct_indic=dim ;
        siz_slct_indic_scenario= dim_scenario;
        dif_incorner_U_L=2;   % decimal acceptance : 2   % delta and pi spliting acceptance : 0.01
        epsilon=.01;
        diffVol_delta=.001;
        diffVol_pi=.001;
        disp(['Num of chosen leftover to solve deter=  ', num2str(dim)]);
        disp(['Num of chosen scenario to solve deter=  ', num2str(dim_scenario)]);
        disp(['Uper-Lower Termination=  ', num2str(epsilon)]);
        disp(['Delta Vol out-in Ter.=   ', num2str(diffVol_delta)]);
        disp(['Pi Vol out-in Ter.=      ', num2str(diffVol_pi)]);
        disp('*Reserve sort= FIFO*');
        disp(['Termination Time=      ', num2str(termint_time),' (s)']);
        disp('--------------------------------');
        %% Basic Inputs
        % Sets
        disp(['Retailers (N)= ', num2str(N)]);
        disp(['Periods (T)=   ', num2str(T)]);
        disp(['Vehicles (K)=  ', num2str(K)]);
        disp(['Breweries (M)= ', num2str(M)]);
        V=M+N;
        % Parameters
        Q= mat.Q;                   % Capacity of vehicle
        FilCap= mat.FilCap;         % Filling capacity
        Cpur= mat.Cpur;             % purchasing cost
        Chyg= mat.Chyg;             % hygienization cost
        Cins= mat.Cins;             % inspection cost
        Cdis= mat.Cdis;             % disposal reimburse
        FixC=mat.FixC;              % fix cost of vehicles
        C= mat.C;                   % variable cost of vehicles
        He= mat.He;
        Hf= mat.Hf;                 % Inventory holding cost of filled bottles
        Lf= mat.Lf;                 % Storage capacity for filled bottles
        Le= mat.Le;                 % Storage capacity for empty bottles
        alpha= mat.alpha;           % Hygienization rate
        Delta(1).A= mat.newA_delta;
        Delta(1).b= mat.newb_delta;
        Pi(1).A= mat.newA_pi;
        Pi(1).b= mat.newb_pi;
        %
        solving_func_time= 0;
        solving_time= 0;
        solving_alg24_func_time= 0;
        solving_Alg24_time= 0;
        deter_time= 0;
        solving_deter_time= 0;
        Modeling_Time=0;
        Alg7_time=0;
        dif_solVSalg24_time=0;
        %
        fuliter=1;
        fuliter_solv_Inner=0;
        fuliter_solv_Outer=0;
        termint_time_indc=0;
        final_Outer_delta=0;
        final_Outer_pi=0;
        indx_final_outer_delta=double.empty(1,0);
        indx_final_outer_pi=double.empty(1,0);
        sol_status_in= 3;
        sol_status_out= 3;
        %
        siz_Total_delta=1;
        siz_Total_pi=1;
        siz_old_acceptbl_delta=1;
        siz_old_acceptbl_pi=1;
        %
        reserv_delta=double.empty(1,0);
        reserv_pi=double.empty(1,0);
        %
        check_Best_diff=Inf;
        %% CLSC Inputs
        % Variables
        purch= sdpvar(M,1);
        m= sdpvar(M,T,'full');
        r= sdpvar(M,T,'full');         % Small "r" for CLSC func
        Iful= sdpvar(M,T+1,'full');    % Iful(:,1) related to time zero
        Iemp= sdpvar(V,T+1,'full');    % Iemp(:,1) related to time zero
        f= sdpvar(N,K,T,'full');
        e= sdpvar(N,K,T,'full');
        v= sdpvar(V,V,K,T,'full');
        u= sdpvar(V,V,K,T,'full');
        X= binvar(V,V,K,T,'full');     % Big "X" for CLSC func
        %
        Iful(1:M,1)= mat.Iful;
        Iemp(1:V,1)= mat.Iemp;
        %
        % For Delta and Delta_scenario cost approximation
        sum_c(1:M,1:K)=0;
        for k=1:K
            for i=1:M
                sum_c(i,k)= sum(C(i,:,k,1),2);
            end
        end
        [val_delta_cost , num_delta_cost]=min(sum_c);
        [~ , indx_K_delta]=min(val_delta_cost);
        indx_M_delta= num_delta_cost(indx_K_delta);
        % For Pi and Pi_scenario cost approximation
        min_c(1:N)=0;
        for j=1:N
            min_c(j)= 2*( min(min(min(C(1:M,j+M,:,:)))) );
            min_He_c(1,T*(j-1)+1:j*T)=min(He(j+M),min_c(j));
        end
        %% constraints
        Constraints = [
            sum(sum(X(1:M,1:M,:,:),3),4) == 0 ;    % Shipment among Breweries is not allowed
            diag(sum(sum(X(M+1:V,M+1:V,:,:),3),4)) == 0 ;    % Routing from a retailer to itself is not allowed
            r <= alpha .* Iemp(1:M,1:T);%Cons 6
            m + r <= repmat(FilCap,1,T); %Cons 7
            Iful(1:M,2:T+1) <= repmat(Lf,1,T);%Cons 9
            Iemp(1:V,2:T+1) <= repmat(Le,1,T);%Cons 10
            sum(X(1:M,M+1:V,1:K,1:T),[1,2]) <= 1;%Cons 11
            sum(sum(X(M+1:V,:,:,1:T),2),3) <= 1;%Cons 12
            sum(X,2) - permute(sum(X,1),[2 1 3 4]) == 0;%Cons 13
            v + u <= Q .* X;%Cons 16
            purch >= 0; m >= 0; r >= 0; Iful(:,2:T+1) >=0;%Cons 17
            Iemp(1:V,2:T+1) >=0;%Cons 18
            f >= 0; e >= 0;%Cons 19
            v >= 0; u >= 0 ];%Cons 20
        for period=1:T
            Constraints = [Constraints;
                Iful(1:M,period) + m(1:M,period) + r(1:M,period) - sum(sum(v(1:M,M+1:V,:,period),2),3) == Iful(1:M,period+1);%Cons 2
                Iemp(1:M,period) + permute(sum(sum(u(M+1:V,1:M,:,period),1),3),[2 1 3 4]) - (r(1:M,period)/alpha)   == Iemp(1:M,period+1); %Cons 3
                m(:,period) <= purch- sum(m(:,1:period-1),2)]; %Cons 8
        end
        for k=1:K
            Constraints = [Constraints;
                permute(sum(v(:,M+1:V,k,1:T),1),[2 1 3 4]) - sum(v(M+1:V,:,k,1:T),2) == f(1:N,k,1:T);%Cons 14
                sum(u(M+1:V,:,k,1:T),2) -  permute(sum(u(:,M+1:V,k,1:T),1),[2 1 3 4]) == e(1:N,k,1:T)];%Cons 15
        end
        % Extra valid constraints
        for period=1:T
            Constraints = [Constraints;
                sum(f(:, :, period), 1)' <= permute(sum( Q(1:M, 1:V, :, period) .* X(1:M, 1:V, :, period), [1 2]),[3,1,2]);   %28_paper
                sum(e(:, :, period), 1)' <= permute(sum( Q(1:M, 1:V, :, period) .* X(1:M, 1:V, :, period), [1 2]),[3,1,2]) ];  %30_paper
            for i=M+1:V
                Constraints = [Constraints;
                    sum(X(i, 1:V, :, period),2) <= sum( X(1:M, 1:V, :, period), [1 2]) ];  %31_paper
            end
        end
        %% Obj Function
        Obj_func= sum(sum(Cpur * purch)) + sum(sum(Chyg * r)) ...
            + sum(sum((Cins-((1-alpha) * Cdis)) * Iemp(1:M,1:T))) ...
            + sum(sum(sum(sum(C .* (v + u))))) ...
            + sum(sum(sum(sum(FixC .* X)))) ...
            + (sum(He' * Iemp(:,2:T+1))) ...
            + (sum(Hf' * Iful(:,2:T+1))) ;
        for tt=1:T
            Obj_func=Obj_func+ He(1:M)'*(purch- sum(m(:,1:tt),2));
        end
        %% Alg 7 Inputs
        prelambda= [];
        lambda(1:2*D,1)=0;
        for i=1:2*D
            if mod(i,2)==0
                lambda(i,1)=Inf;
            else
                lambda(i,1)=-Inf;
            end
        end
        %% Alg 24 Inputs
        x=sdpvar(D,1);  % Small "x" for Alg24 func
        PoolscenarioU_delta=[];
        PoolscenarioU_pi=[];
        %% End of Main Inputs
        %% Pre iteration
        disp('#################################');
        disp(['Iteration (', num2str(fuliter),') ...'])
        disp('==============(Delta)==============');
        disp('finding box');
        iter_solv_Inner=0;
        iter_solv_Outer=0;
        Objective(1).inner(1,1)= -Inf;
        Objective(1).inner(1,2:5)= -1;
        Objective(1).outer(1,1)= Inf;
        Objective(1).outer(1,2:3)= -1;
        teek0=tic;
        teek=tic;
        [ Delta(1).outcornerU(1:D,1), Delta(1).incornerL(1:D,1), Delta(1).incornerU(1:D,1), Delta(1).scenario(1:D,:), PoolscenarioU_delta(1:D,:), solving_alg24_func_time, solving_func_time, Delta(1).Volin, Delta(1).Volout] = Alg24_func (Delta(1).A, Delta(1).b, x, D, PoolscenarioU_delta); % Alg 24
        Modeling_Time=Modeling_Time+toc(teek);
        Delta(1).forbiden_pi=[];
        solving_Alg24_time= solving_Alg24_time+solving_alg24_func_time;
        solving_alg24_func_time= 0;
        solving_time= solving_time + solving_func_time;
        solving_func_time= 0;
        siz_scenarioU_delta(1,1)=size(Delta(1).scenario,2);
        dif_out_in_delta(1,1)= Delta(1).Volout - Delta(1).Volin;
        % first delta selection
        siz_slct_delta(fuliter,1)=1;                                                      % #size delta select (solve)
        siz_slct_delta_scenario(fuliter,1)= min(siz_slct_indic_scenario, siz_scenarioU_delta(1,1));               % #size delta scenario select (solve)
        Delta(fuliter).chosen_calc= ones(1,siz_slct_delta_scenario(fuliter,1)+1);  % delta and delta scenarios
        Delta(fuliter+1).chosen_split= 1;
        cost_delta_scenario(1, 1:D)= 0;
        if siz_scenarioU_delta(1,1) ~= 0
            ind_row=1;
            for ind_N=M+1:V
                cost_delta_scenario(1, 1: siz_scenarioU_delta(1,1))= cost_delta_scenario(1, 1: siz_scenarioU_delta(1,1))+ sum((2* C(indx_M_delta,ind_N,indx_K_delta,1)* Delta(1).scenario(ind_row: ind_row +T-1, 1: siz_scenarioU_delta(1,1))),1);
                ind_row= ind_row+T;
            end
            [~, indx_delta_scenario]= sort(cost_delta_scenario,'descend');
            Delta(fuliter).chosen_scenario_colmn= [0, indx_delta_scenario(1:siz_slct_delta_scenario(fuliter,1))];  % column(s) of delta scenario
        end
        disp('================(Pi)===============');
        disp('finding box');
        teek=tic;
        [ Pi(1).outcornerU(1:D,1), Pi(1).incornerL(1:D,1), Pi(1).incornerU(1:D,1), Pi(1).scenario(1:D,:), PoolscenarioU_pi(1:D,:), solving_alg24_func_time, solving_func_time, Pi(1).Volin, Pi(1).Volout] = Alg24_func (Pi(1).A, Pi(1).b, x, D, PoolscenarioU_pi); % Alg 24
        Modeling_Time=Modeling_Time+toc(teek);
        Pi(1).tree=1;
        solving_Alg24_time= solving_Alg24_time+solving_alg24_func_time;
        solving_alg24_func_time= 0;
        solving_time= solving_time + solving_func_time;
        solving_func_time= 0;
        siz_scenarioU_pi(1,1)=size( Pi(1).scenario,2);
        dif_out_in_pi(1,1)= Pi(1).Volout- Pi(1).Volin;
        % first pi selection
        siz_slct_pi(fuliter,1)=1;                                                        % #size pi select (solve)
        siz_slct_pi_scenario(fuliter,1)= min(siz_slct_indic_scenario, siz_scenarioU_pi(1,1));                    % #size pi scenario select (solve)
        Pi(fuliter).chosen_calc= ones(1,siz_slct_pi_scenario(fuliter,1)+1);  % pi and pi scenarios
        Pi(fuliter+1).chosen_split= 1;
        cost_pi_scenario(1, 1:D)=0;
        if siz_scenarioU_pi(1,1) ~= 0
            cost_pi_scenario(1, 1: siz_scenarioU_pi(1,1))= min_He_c * Pi(1).scenario;
            [~, indx_pi_scenario]= sort(cost_pi_scenario,'descend');
            Pi(fuliter).chosen_scenario_colmn= [0, indx_pi_scenario(1:siz_slct_pi_scenario(fuliter,1))];           % column(s) of pi scenario
        end
        disp('*********************************');
        %% Pre iteration solving
        disp('(Inner, Outer) new delta <-> new pi ...');
        Best.incornerU_delta= Inf;
        Best.incornerU_pi= Inf;
        for i=1:siz_slct_delta_scenario(fuliter,1)+1
            for j=1:siz_slct_pi_scenario(fuliter,1)+1
                %
                fultime= toc(teek0);
                solution_time= fultime-Modeling_Time+solving_time;
                termint_time_func= max((termint_time-solution_time),0);
                %
                iter_solv_Inner= iter_solv_Inner+1;
                teek=tic;
                if i == 1
                    if j == 1 % both inners
                        [Objective(fuliter).inner(iter_solv_Inner,1), deter_time, solving_func_time, sol_status_in, Objective(fuliter).inner_var(iter_solv_Inner)] = CLSC_func (Delta(1).incornerU, Pi(1).incornerU, Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
                        Objective(fuliter).inner(iter_solv_Inner,2:5)=[1, 1, 0, 0];
                    else % delta inners, pi scenarios
                        [Objective(fuliter).inner(iter_solv_Inner,1), deter_time, solving_func_time, sol_status_in, Objective(fuliter).inner_var(iter_solv_Inner)] = CLSC_func (Delta(1).incornerU, Pi(1).scenario(:,Pi(fuliter).chosen_scenario_colmn(j)), Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
                        Objective(fuliter).inner(iter_solv_Inner,2:5)=[1, 1, 0, Pi(fuliter).chosen_scenario_colmn(j)];    % [ delta, pi, scenario_delta, scenario_pi ]
                    end
                else
                    if j == 1 % delta scenarios, pi inners
                        [Objective(fuliter).inner(iter_solv_Inner,1), deter_time, solving_func_time, sol_status_in, Objective(fuliter).inner_var(iter_solv_Inner)] = CLSC_func (Delta(1).scenario(:,Delta(fuliter).chosen_scenario_colmn(i)), Pi(1).incornerU, Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
                        Objective(fuliter).inner(iter_solv_Inner,2:5)=[1, 1, Delta(fuliter).chosen_scenario_colmn(i), 0];
                    else % both scenarios
                        [Objective(fuliter).inner(iter_solv_Inner,1), deter_time, solving_func_time, sol_status_in, Objective(fuliter).inner_var(iter_solv_Inner)] = CLSC_func (Delta(1).scenario(:,Delta(fuliter).chosen_scenario_colmn(i)), Pi(1).scenario(:,Pi(fuliter).chosen_scenario_colmn(j)), Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
                        Objective(fuliter).inner(iter_solv_Inner,2:5)=[1, 1, Delta(fuliter).chosen_scenario_colmn(i), Pi(fuliter).chosen_scenario_colmn(j)];    % [ delta, pi, scenario_delta, scenario_pi ]
                    end
                end
                Modeling_Time=Modeling_Time+toc(teek);
                solving_deter_time= solving_deter_time +deter_time;
                deter_time= 0;
                solving_time= solving_time + solving_func_time;
                solving_func_time= 0;
                if sol_status_in ~= 0        % it could not be solved because of time termination
                    disp(['sol_status_in= ', num2str(sol_status_in)]);
                    disp(['Solving Time= ', num2str(solving_time)]);
                    disp(['Solution Time= ', num2str(solution_time)]);
                    disp(['Deterministic Time= ', num2str(solving_deter_time)]);
                    Objective(fuliter).inner(iter_solv_Inner,1)= -Inf;
                    Objective(fuliter).inner(iter_solv_Inner,2:5)= -1;
                    iter_solv_Inner= iter_solv_Inner-1;
                    termint_time_indc=2;
                    break
                end
            end
            if termint_time_indc==2
                break
            end
        end
        %
        fultime= toc(teek0);
        solution_time= fultime-Modeling_Time+solving_time;
        termint_time_func= max((termint_time-solution_time),0);
        %
        Best.outcornerU_delta= Inf;
        Best.outcornerU_pi= Inf;
        if termint_time_indc~=2
            iter_solv_Outer= iter_solv_Outer+1;
            teek=tic;
            [Objective(fuliter).outer(iter_solv_Outer,1), deter_time, solving_func_time, sol_status_out, Objective(fuliter).outer_var(iter_solv_Outer)] = CLSC_func (Delta(1).outcornerU, Pi(1).outcornerU, Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
            Modeling_Time=Modeling_Time+toc(teek);
            solving_deter_time= solving_deter_time +deter_time;
            deter_time= 0;
            solving_time= solving_time + solving_func_time;
            solving_func_time= 0;
            if sol_status_out == 0
                Objective(fuliter).outer(iter_solv_Outer,2:3)=[1, 1];
                Best.outcornerU_delta= Delta(1).outcornerU;
                Best.outcornerU_pi= Pi(1).outcornerU;
            else        % it could not be solved because of time termination
                disp(['sol_status_out= ', num2str(sol_status_out)])
                disp(['Solving Time= ', num2str(solving_time)]);
                disp(['Solution Time= ', num2str(solution_time)]);
                disp(['Deterministic Time= ', num2str(solving_deter_time)]);
                Objective(fuliter).outer(iter_solv_Outer,1)= Inf;
                Objective(fuliter).outer(iter_solv_Outer,2:3)= -1;
                iter_solv_Outer= iter_solv_Outer-1;
                termint_time_indc=2;
            end
        end
        [maxin(fuliter,1), rowmaxin]= max(Objective(fuliter).inner(:,1));
        maxout(fuliter,1)= Objective(fuliter).outer(1,1);
        out_in_differnc(fuliter)= maxout(fuliter,1) - maxin(fuliter,1);
        Best.Differenc=out_in_differnc(1);
        best_clmn_indc_delta_scenario=Objective(fuliter).inner(rowmaxin,4);
        best_clmn_indc_pi_scenario=Objective(fuliter).inner(rowmaxin,5);
        Best.Inner_iter_row=[fuliter, rowmaxin];
        Best.Outer_iter_row=[fuliter, 1];
        disp('----------------------------------')
        if best_clmn_indc_delta_scenario == 0
            disp('Best Inner: delta [1]')
            Best.incornerU_delta= Delta(1).incornerU;
            Delta(fuliter).chosen_calc= 1;
            siz_slct_delta_scenario(fuliter,1) =0;
            Delta(fuliter).chosen_scenario_colmn =0;
        elseif best_clmn_indc_delta_scenario ~= -1
            disp(['Best Inner: delta [1] scenario(',num2str(best_clmn_indc_delta_scenario),')'])
            Best.incornerU_delta= Delta(1).scenario(:,best_clmn_indc_delta_scenario);
            siz_slct_delta(fuliter,1) =0;
            Delta(fuliter).chosen_calc =1;
            siz_slct_delta_scenario(fuliter,1) =1;
            Delta(fuliter).chosen_scenario_colmn = best_clmn_indc_delta_scenario;
        end
        if best_clmn_indc_pi_scenario == 0
            disp('Best Inner: pi [1]')
            Best.incornerU_pi= Pi(1).incornerU;
            Pi(fuliter).chosen_calc= 1;
            siz_slct_pi_scenario(fuliter,1)= 0;
            Pi(fuliter).chosen_scenario_colmn= 0;
        elseif best_clmn_indc_pi_scenario ~= -1
            disp(['Best Inner: pi [1] scenario(',num2str(best_clmn_indc_pi_scenario),')'])
            Best.incornerU_pi= Pi(1).scenario(:,best_clmn_indc_pi_scenario);
            siz_slct_pi(fuliter,1)= 0;
            Pi(fuliter).chosen_calc= 1;
            siz_slct_pi_scenario(fuliter,1)= 1;
            Pi(fuliter).chosen_scenario_colmn= best_clmn_indc_pi_scenario;
        end
        Best.Inner= maxin(fuliter,1);
        Best.Outer= maxout(fuliter,1);
        fuliter_solv_Inner= fuliter_solv_Inner +iter_solv_Inner;
        fuliter_solv_Outer= fuliter_solv_Outer +iter_solv_Outer;
        disp('----------------------------------')
        disp(['Total maxin= ', num2str(Best.Inner)]);
        disp(['Row Maxin=  ',num2str(rowmaxin)]);
        disp(['Total maxout= ', num2str(Best.Outer)]);
        disp(['Row Maxout=  ',num2str(1)]);
        disp(['Min difference= ', num2str(Best.Differenc)]);
        disp('----------------------------------');
        fultime= toc(teek0);
        solution_time= fultime-Modeling_Time+solving_time;
        disp(['Full Time=  ', num2str(fultime)]);
        disp(['Modeling Time= ', num2str(Modeling_Time)]);
        disp(['Solution Time= ', num2str(solution_time)]);
        disp(['Deterministic Time= ', num2str(solving_deter_time)]);
        disp(['Alg24 solver Time= ', num2str(solving_Alg24_time)]);
        disp(['Number of inner solving: ', num2str(fuliter_solv_Inner)]);
        disp(['Number of outer solving: ', num2str(fuliter_solv_Outer)]);
        disp(['Number of All solving: ', num2str(fuliter_solv_Inner+fuliter_solv_Outer)]);
        if termint_time_indc==2
            disp('*********************************')
            disp('Termination: Time termination')
        end
        %% End of Pre iteration
        %% Main loop
        while out_in_differnc(fuliter) > epsilon && termint_time_indc==0
            fuliter=fuliter+1;
            disp('#################################');
            disp(['Iteration (', num2str(fuliter),') ...']);
            %% Set for delta and pi
            Delta(fuliter).chosen_calc= double.empty(1,0);
            Pi(fuliter).chosen_calc= double.empty(1,0);
            siz_acceptbl_delta= 0;
            siz_acceptbl_pi= 0;
            cost_delta= zeros(1);                        % all previous cost_delta will be zero
            cost_delta_scenario= zeros(1);
            cost_pi= zeros(1);
            cost_pi_scenario= zeros(1);
            siz_slct_delta(fuliter,1)= 0;
            siz_slct_pi(fuliter,1)= 0;
            siz_slct_delta_scenario(fuliter,1)=0;
            siz_slct_pi_scenario(fuliter,1)=0;
            indcat_split_delta=0;
            indcat_split_pi=0;
            %% Delta iteration
            disp('==============(Delta)==============');
            if siz_old_acceptbl_delta~=0
                disp(['siz_Total_delta until now= {', num2str(siz_Total_delta),'}']);
                disp('- - - - - - - - - - -');
                while indcat_split_delta == 0
                    iter_split_delta= Delta(fuliter).chosen_split;
                    if  termint_time_indc==2
                        break
                    end
                    if dif_out_in_delta(iter_split_delta,1) > diffVol_delta
                        teek=tic;
                        [newA, newb, solving_func_time]= Alg7_func (Delta(iter_split_delta).A, Delta(iter_split_delta).b, D, Delta(iter_split_delta).incornerL, Delta(iter_split_delta).incornerU, prelambda, lambda); % Alg 7
                        Modeling_Time=Modeling_Time+toc(teek);
                        solving_time= solving_time + solving_func_time;
                        Alg7_time=Alg7_time+solving_func_time;
                        %                         disp(['Alg7 time: ',num2str(Alg7_time)]);
                        solving_func_time= 0;
                        for i= 1:size(newA,2)
                            fultime= toc(teek0);
                            solution_time= fultime-Modeling_Time+solving_time;
                            if solution_time>= termint_time
                                termint_time_indc=2;
                                disp('Termination: Time termination');
                                break
                            end
                            siz_Total_delta= siz_Total_delta+ 1;
                            teek=tic;
                            [ Delta(siz_Total_delta).outcornerU, Delta(siz_Total_delta).incornerL, Delta(siz_Total_delta).incornerU, Delta(siz_Total_delta).scenario, PoolscenarioU_delta, solving_alg24_func_time, solving_func_time, Delta(siz_Total_delta).Volin, Delta(siz_Total_delta).Volout] = Alg24_func (newA(i).part, newb(i).part, x, D, PoolscenarioU_delta); % Alg 24
                            Modeling_Time=Modeling_Time+toc(teek);
                            solving_Alg24_time= solving_Alg24_time+solving_alg24_func_time;
                            solving_time= solving_time + solving_func_time;
                            dif_solVSalg24_time=dif_solVSalg24_time+(solving_func_time-solving_alg24_func_time);
                            solving_alg24_func_time= 0;
                            solving_func_time= 0;
                            if ismember(0, round(Delta(siz_Total_delta).incornerU - Delta(siz_Total_delta).incornerL,dif_incorner_U_L)) == 0
                                siz_scenarioU_delta(siz_Total_delta,1)=size( Delta(siz_Total_delta).scenario,2);
                                siz_acceptbl_delta= siz_acceptbl_delta+ 1;
                                Delta(siz_Total_delta).forbiden_pi=Delta(iter_split_delta).forbiden_pi;
                                %
                                Delta(siz_Total_delta).A= newA(i).part;
                                Delta(siz_Total_delta).b= newb(i).part;
                                dif_out_in_delta(siz_Total_delta,1)= Delta(siz_Total_delta).Volout-Delta(siz_Total_delta).Volin; % it's on Totalsize delta
                                %
                                cost_delta(1, siz_acceptbl_delta)=0; % cost delta calculation
                                ind_row=1;
                                for ind_N=M+1:V
                                    cost_delta(1, siz_acceptbl_delta)= cost_delta(1, siz_acceptbl_delta)+ sum(2* C(indx_M_delta,ind_N,indx_K_delta,1)* Delta(siz_Total_delta).incornerU(ind_row: ind_row +T-1, 1));
                                    ind_row= ind_row+T;
                                end
                                if siz_scenarioU_delta(siz_Total_delta,1) ~=0  % cost delta scenario calculation
                                    cost_delta_scenario(siz_acceptbl_delta, 1:siz_scenarioU_delta(siz_Total_delta,1))= 0;
                                    ind_row=1;
                                    for ind_N=M+1:V
                                        cost_delta_scenario(siz_acceptbl_delta, 1:siz_scenarioU_delta(siz_Total_delta,1))= cost_delta_scenario(siz_acceptbl_delta, 1:siz_scenarioU_delta(siz_Total_delta,1))+ sum((2* C(indx_M_delta,ind_N,indx_K_delta,1)* Delta(siz_Total_delta).scenario(ind_row: ind_row +T-1, 1:siz_scenarioU_delta(siz_Total_delta,1))),1);
                                        ind_row= ind_row+T;
                                    end
                                end
                            else
                                siz_Total_delta= siz_Total_delta- 1;
                                disp('"One new Delat polytope rejected"');
                            end
                        end
                        disp(['Delta(', num2str(iter_split_delta),'): "', num2str(dif_out_in_delta(iter_split_delta,1)),'", +',num2str(siz_acceptbl_delta),' = ', num2str(siz_Total_delta)]);
                    else
                        final_Outer_delta= final_Outer_delta+1;
                        indx_final_outer_delta(1,final_Outer_delta)= iter_split_delta;  % for OuterU Delta belong to small polyttops
                        disp(['Delta(', num2str(iter_split_delta),'): "', num2str(dif_out_in_delta(iter_split_delta,1)),'" -> final_Outer_delta (',num2str(final_Outer_delta),')']);
                    end
                    if siz_acceptbl_delta == 0     % no new polytope because of shape
                        if size(reserv_delta,2) ~=0
                            disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
                            Delta(fuliter).chosen_split= reserv_delta(1); % FIFO
                            disp(['chosen split delta changed to Reserve: [', num2str(reserv_delta(1)) ,']']);
                            reserv_delta(1)=[];
                        else
                            indcat_split_delta=1;
                            %                             if final_Outer_delta==0
                            %                                 final_Outer_delta= final_Outer_delta+1;
                            %                                 indx_final_outer_delta(1,final_Outer_delta)= iter_split_delta;  % for OuterU Delta belong to small polyttops
                            %                                 disp(['Delta(', num2str(iter_split_delta),'): "', num2str(dif_out_in_delta(iter_split_delta,1)),'" -> final_Outer_delta (',num2str(final_Outer_delta),')']);
                            %                             end
                        end
                    else
                        indcat_split_delta=1;
                    end
                end
                if  termint_time_indc==2
                    break
                end
                siz_old_acceptbl_delta= siz_acceptbl_delta; % updating sizeoldA_delta for next Step
                if siz_old_acceptbl_delta ~= 0
                    siz_slct_delta(fuliter,1)=min(siz_slct_indic, siz_old_acceptbl_delta);        % #size delta select (solve)
                    [~, indx_cost_delta]= sort(cost_delta,'descend');    % to sort Delta polytops by cost for next step
                    Delta(fuliter).chosen_calc(1, 1:siz_slct_delta(fuliter,1))= siz_Total_delta-siz_old_acceptbl_delta+indx_cost_delta(1:siz_slct_delta(fuliter,1));  % deltas wil be selected (without siz_choos_byCost_delta)
                    Delta(fuliter).chosen_scenario_colmn(1, 1:siz_slct_delta(fuliter,1))= zeros(1,siz_slct_delta(fuliter,1));
                    disp('- - - - - - - - - - -');
                    % delta_scenario_choose
                    siz_slct_delta_scenario(fuliter,1)=min(siz_slct_indic_scenario, sum(siz_scenarioU_delta(siz_Total_delta-siz_old_acceptbl_delta+1:siz_Total_delta)));         % #size delta scenario select (solve)
                    [~, indx_delta_scenario]= sort(cost_delta_scenario(:),'descend');
                    divd_val= indx_delta_scenario(1:siz_slct_delta_scenario(fuliter,1))'/(siz_acceptbl_delta);
                    Delta(fuliter).chosen_scenario_colmn(1, end+1:end+siz_slct_delta_scenario(fuliter,1))= ceil(divd_val);                          % column(s) of delta scenario
                    decmal =divd_val- fix(divd_val);
                    Delta(fuliter).chosen_calc(1, end+1:end+siz_slct_delta_scenario(fuliter,1))= round(decmal*siz_acceptbl_delta);
                    Delta(fuliter).chosen_calc(1, Delta(fuliter).chosen_calc(1,:)==0)= siz_acceptbl_delta;
                    Delta(fuliter).chosen_calc(1, siz_slct_delta(fuliter,1)+1:end)= siz_Total_delta-siz_old_acceptbl_delta+Delta(fuliter).chosen_calc(1, siz_slct_delta(fuliter,1)+1:end);
                    disp(['delta to calc    : [', num2str(Delta(fuliter).chosen_calc) ,']'])
                    disp(['scenario to calc : [', num2str(Delta(fuliter).chosen_scenario_colmn) ,']'])
                else
                    disp('* No new Delta from split');
                end
            else
                disp('* No new Delta from split');
            end
            %% End of Delta iteration
            %% pi iteration
            disp('================(Pi)===============');
            if siz_old_acceptbl_pi~=0
                disp(['siz_Total_pi until now= {', num2str(siz_Total_pi),'}']);
                disp('- - - - - - - - - - -');
                while indcat_split_pi == 0
                    iter_split_pi= Pi(fuliter).chosen_split;
                    if  termint_time_indc==2
                        break
                    end
                    if dif_out_in_pi(iter_split_pi,1) > diffVol_pi
                        teek=tic;
                        [newA, newb, solving_func_time]= Alg7_func (Pi(iter_split_pi).A, Pi(iter_split_pi).b, D, Pi(iter_split_pi).incornerL, Pi(iter_split_pi).incornerU, prelambda, lambda); % Alg 7
                        Modeling_Time=Modeling_Time+toc(teek);
                        solving_time= solving_time + solving_func_time;
                        Alg7_time=Alg7_time+solving_func_time;
                        %                         disp(['Alg7 time: ',num2str(Alg7_time)]);
                        solving_func_time= 0;
                        for i= 1:size(newA,2)
                            fultime= toc(teek0);
                            solution_time= fultime-Modeling_Time+solving_time;
                            if solution_time>= termint_time
                                termint_time_indc=2;
                                disp('Termination: Time termination');
                                break
                            end
                            siz_Total_pi= siz_Total_pi+ 1;
                            teek=tic;
                            [ Pi(siz_Total_pi).outcornerU, Pi(siz_Total_pi).incornerL, Pi(siz_Total_pi).incornerU, Pi(siz_Total_pi).scenario, PoolscenarioU_pi, solving_alg24_func_time, solving_func_time, Pi(siz_Total_pi).Volin, Pi(siz_Total_pi).Volout] = Alg24_func (newA(i).part, newb(i).part, x, D, PoolscenarioU_pi); % Alg 24
                            Modeling_Time=Modeling_Time+toc(teek);
                            solving_Alg24_time= solving_Alg24_time+solving_alg24_func_time;
                            solving_time= solving_time + solving_func_time;
                            dif_solVSalg24_time=dif_solVSalg24_time+(solving_func_time-solving_alg24_func_time);
                            solving_alg24_func_time= 0;
                            solving_func_time= 0;
                            if ismember(0, round(Pi(siz_Total_pi).incornerU - Pi(siz_Total_pi).incornerL,dif_incorner_U_L)) == 0
                                siz_scenarioU_pi(siz_Total_pi,1)=size( Pi(siz_Total_pi).scenario,2);
                                siz_acceptbl_pi= siz_acceptbl_pi+ 1;
                                Pi(siz_Total_pi).tree=[Pi(iter_split_pi).tree, siz_Total_pi];
                                %
                                Pi(siz_Total_pi).A= newA(i).part;
                                Pi(siz_Total_pi).b= newb(i).part;
                                dif_out_in_pi(siz_Total_pi,1)= Pi(siz_Total_pi).Volout-Pi(siz_Total_pi).Volin;  % it's on Totalsize pi
                                %
                                cost_pi(1,siz_acceptbl_pi)= min_He_c * Pi(siz_Total_pi).incornerU;
                                if siz_scenarioU_pi(siz_Total_pi,1) ~= 0
                                    cost_pi_scenario(siz_acceptbl_pi, 1:siz_scenarioU_pi(siz_Total_pi,1))= min_He_c * Pi(siz_Total_pi).scenario;
                                end
                            else
                                siz_Total_pi= siz_Total_pi- 1;
                                disp('"One new Pi polytope rejected"');
                            end
                        end
                        disp(['Pi(', num2str(iter_split_pi),'): "', num2str(dif_out_in_pi(iter_split_pi,1)),'", +',num2str(siz_acceptbl_pi),' = ', num2str(siz_Total_pi)]);
                    else
                        final_Outer_pi= final_Outer_pi+1;
                        indx_final_outer_pi(1,final_Outer_pi)= iter_split_pi;  % for OuterU pi belong to small polyttops
                        disp(['Pi(', num2str(iter_split_pi),'): "', num2str(dif_out_in_pi(iter_split_pi,1)),'" -> final_Outer_pi (',num2str(final_Outer_pi),')']);
                    end
                    if siz_acceptbl_pi == 0
                        if size(reserv_pi,2) ~=0
                            disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
                            Pi(fuliter).chosen_split= reserv_pi(1); % FIFO
                            disp(['chosen split pi changed to Reserve: [', num2str(reserv_pi(1)) ,']']);
                            reserv_pi(1)=[];
                        else
                            indcat_split_pi=1;
                            %                             if final_Outer_pi==0
                            %                                 final_Outer_pi= final_Outer_pi+1;
                            %                                 indx_final_outer_pi(1,final_Outer_pi)= iter_split_pi;  % for OuterU pi belong to small polyttops
                            %                                 disp(['Pi(', num2str(iter_split_pi),'): "', num2str(dif_out_in_pi(iter_split_pi,1)),'" -> final_Outer_pi (',num2str(final_Outer_pi),')']);
                            %                             end
                        end
                    else
                        indcat_split_pi=1;
                    end
                end
                if  termint_time_indc==2
                    break
                end
                siz_old_acceptbl_pi= siz_acceptbl_pi; % updating sizeoldA_pi for next Step
                if siz_old_acceptbl_pi ~= 0
                    siz_slct_pi(fuliter,1)=min(siz_slct_indic, siz_old_acceptbl_pi);        % #size pi select (solve)
                    [~, indx_cost_pi]= sort(cost_pi,'descend');    % to sort pi polytops by cost for next step
                    Pi(fuliter).chosen_calc(1, 1:siz_slct_pi(fuliter,1))= siz_Total_pi-siz_old_acceptbl_pi+indx_cost_pi(1:siz_slct_pi(fuliter,1));  % pis wil be selected (without siz_choos_byCost_pi)
                    Pi(fuliter).chosen_scenario_colmn(1, 1:siz_slct_pi(fuliter,1))= zeros(1,siz_slct_pi(fuliter,1));
                    disp('- - - - - - - - - - -');
                    % pi_scenario_choos
                    siz_slct_pi_scenario(fuliter,1)=min(siz_slct_indic_scenario, sum(siz_scenarioU_pi(siz_Total_pi-siz_old_acceptbl_pi+1:siz_Total_pi)));         % #size pi scenario select (solve)
                    [~, indx_pi_scenario]= sort(cost_pi_scenario(:),'descend');
                    divd_val= indx_pi_scenario(1:siz_slct_pi_scenario(fuliter,1))'/(siz_acceptbl_pi);
                    Pi(fuliter).chosen_scenario_colmn(1, end+1:end+siz_slct_pi_scenario(fuliter,1))= ceil(divd_val);                          % column(s) of pi scenario
                    decmal =divd_val- fix(divd_val);
                    Pi(fuliter).chosen_calc(1, end+1:end+siz_slct_pi_scenario(fuliter,1))= round(decmal*siz_acceptbl_pi);
                    Pi(fuliter).chosen_calc(1, Pi(fuliter).chosen_calc(1,:)==0)= siz_acceptbl_pi;
                    Pi(fuliter).chosen_calc(1, siz_slct_pi(fuliter,1)+1:end)= siz_Total_pi-siz_old_acceptbl_pi+Pi(fuliter).chosen_calc(1, siz_slct_pi(fuliter,1)+1:end);
                    disp(['pi to calc       : [', num2str(Pi(fuliter).chosen_calc) ,']'])
                    disp(['scenario to calc : [', num2str(Pi(fuliter).chosen_scenario_colmn) ,']'])
                else
                    disp('* No new Pi from split');
                end
            else
                disp('* No new Pi from split');
            end
            %% End of pi iteration
            %% Difference and calculation
            Delta(fuliter+1).chosen_split=double.empty(1,0);
            Pi(fuliter+1).chosen_split=double.empty(1,0);
            if siz_old_acceptbl_delta ~=0 || siz_old_acceptbl_pi ~=0
                %% Calculate Obj func
                disp('*********************************');
                iter_solv_Inner=0;
                iter_solv_Outer=0;
                Objective(fuliter).inner(1,1)= -Inf;
                Objective(fuliter).inner(1,2:5)= -1;
                Objective(fuliter).outer(1,1)= -Inf;
                Objective(fuliter).outer(1,2:3)= -1;
                for calc_step= 1:3
                    if  termint_time_indc==1
                        break
                    end
                    indcat_calc_step=0;
                    if calc_step == 1 && siz_old_acceptbl_delta ~=0 && siz_old_acceptbl_pi ~=0 % (new delta , new pi)
                        indx_calc_delta_inner=fuliter;
                        indx_calc_pi_inner=fuliter;
                        indx_calc_delta_outer=Delta(fuliter).chosen_calc(1:siz_slct_delta(fuliter,1));
                        indx_calc_pi_outer=Pi(fuliter).chosen_calc(1:siz_slct_pi(fuliter,1));
                        indcat_calc_step=1;
                        disp('(Inner, Outer) new delta <-> new pi ...')
                    elseif calc_step == 2 && siz_old_acceptbl_delta ~=0 % (new delta , old pi)
                        indx_calc_delta_inner=fuliter;
                        indx_calc_pi_inner=1:fuliter-1;
                        indx_calc_delta_outer=Delta(fuliter).chosen_calc(1:siz_slct_delta(fuliter,1));
                        indx_calc_pi_outer=indx_final_outer_pi;
                        indcat_calc_step=1;
                        disp('(Inner, Outer) new delta <-> old pi ...')
                    elseif calc_step == 3 && siz_old_acceptbl_pi ~=0 % (old delta , new pi)
                        indx_calc_delta_inner=1:fuliter-1;
                        indx_calc_pi_inner=fuliter;
                        indx_calc_delta_outer=indx_final_outer_delta;
                        indx_calc_pi_outer=Pi(fuliter).chosen_calc(1:siz_slct_pi(fuliter,1));
                        indcat_calc_step=1;
                        disp('(Inner, Outer) old delta <-> new pi ...')
                    end
                    if indcat_calc_step ==1
                        for k_delta=indx_calc_delta_inner
                            if  termint_time_indc==1
                                break
                            end
                            for i= 1:siz_slct_delta(k_delta,1)+siz_slct_delta_scenario(k_delta,1)
                                if  termint_time_indc==1
                                    break
                                end
                                for k_pi=indx_calc_pi_inner
                                    if  termint_time_indc==1
                                        break
                                    end
                                    for j= 1:siz_slct_pi(k_pi,1)+siz_slct_pi_scenario(k_pi,1)
                                        fultime= toc(teek0);
                                        solution_time= fultime-Modeling_Time+solving_time;
                                        if solution_time>= termint_time
                                            termint_time_indc=1;
                                            break
                                        end
                                        %
                                        indcat_forbiden=0;
                                        if size(unique([Pi(Pi(k_pi).chosen_calc(j)).tree, Delta(Delta(k_delta).chosen_calc(i)).forbiden_pi]),2) < size([Pi(Pi(k_pi).chosen_calc(j)).tree, Delta(Delta(k_delta).chosen_calc(i)).forbiden_pi],2)
                                            indcat_forbiden=1;
                                        end
                                        %
                                        if indcat_forbiden == 0
                                            fultime= toc(teek0);
                                            solution_time= fultime-Modeling_Time+solving_time;
                                            if solution_time>= termint_time
                                                termint_time_indc=1;
                                                break
                                            else
                                                termint_time_func= termint_time-solution_time;
                                            end
                                            iter_solv_Inner= iter_solv_Inner+1;
                                            teek=tic;
                                            if i<=siz_slct_delta(k_delta,1)
                                                if j<=siz_slct_pi(k_pi,1)
                                                    [Objective(fuliter).inner(iter_solv_Inner,1), deter_time, solving_func_time, sol_status_in, Objective(fuliter).inner_var(iter_solv_Inner)] = CLSC_func (Delta(Delta(k_delta).chosen_calc(i)).incornerU, Pi(Pi(k_pi).chosen_calc(j)).incornerU, Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
                                                    Objective(fuliter).inner(iter_solv_Inner,2:5)=[Delta(k_delta).chosen_calc(i), Pi(k_pi).chosen_calc(j), 0, 0];
                                                else
                                                    [Objective(fuliter).inner(iter_solv_Inner,1), deter_time, solving_func_time, sol_status_in, Objective(fuliter).inner_var(iter_solv_Inner)] = CLSC_func (Delta(Delta(k_delta).chosen_calc(i)).incornerU, Pi(Pi(k_pi).chosen_calc(j)).scenario(:,Pi(k_pi).chosen_scenario_colmn(j)), Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
                                                    Objective(fuliter).inner(iter_solv_Inner,2:5)=[Delta(k_delta).chosen_calc(i), Pi(k_pi).chosen_calc(j), 0, Pi(k_pi).chosen_scenario_colmn(j)];
                                                end
                                            else
                                                if j<=siz_slct_pi(k_pi,1)
                                                    [Objective(fuliter).inner(iter_solv_Inner,1), deter_time, solving_func_time, sol_status_in, Objective(fuliter).inner_var(iter_solv_Inner)] = CLSC_func (Delta(Delta(k_delta).chosen_calc(i)).scenario(:,Delta(k_delta).chosen_scenario_colmn(i)), Pi(Pi(k_pi).chosen_calc(j)).incornerU, Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
                                                    Objective(fuliter).inner(iter_solv_Inner,2:5)=[Delta(k_delta).chosen_calc(i), Pi(k_pi).chosen_calc(j), Delta(k_delta).chosen_scenario_colmn(i), 0];
                                                else
                                                    [Objective(fuliter).inner(iter_solv_Inner,1), deter_time, solving_func_time, sol_status_in, Objective(fuliter).inner_var(iter_solv_Inner)] = CLSC_func (Delta(Delta(k_delta).chosen_calc(i)).scenario(:,Delta(k_delta).chosen_scenario_colmn(i)), Pi(Pi(k_pi).chosen_calc(j)).scenario(:,Pi(k_pi).chosen_scenario_colmn(j)), Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
                                                    Objective(fuliter).inner(iter_solv_Inner,2:5)=[Delta(k_delta).chosen_calc(i), Pi(k_pi).chosen_calc(j), Delta(k_delta).chosen_scenario_colmn(i), Pi(k_pi).chosen_scenario_colmn(j)];
                                                end
                                            end
                                            Modeling_Time=Modeling_Time+toc(teek);
                                            solving_deter_time= solving_deter_time +deter_time;
                                            deter_time= 0;
                                            solving_time= solving_time + solving_func_time;
                                            solving_func_time= 0;
                                            if sol_status_in ~= 0
                                                disp(['sol_status_in= ', num2str(sol_status_in)]);
                                                disp(['Solving Time= ', num2str(solving_time)]);
                                                disp(['Solution Time= ', num2str(solution_time)]);
                                                disp(['Deterministic Time= ', num2str(solving_deter_time)]);
                                                fultime= toc(teek0);
                                                disp(['Full Time=  ', num2str(fultime)]);
                                                disp(['Modeling Time= ', num2str(Modeling_Time)]);
                                                Objective(fuliter).inner(iter_solv_Inner,1)= -Inf;
                                                Objective(fuliter).inner(iter_solv_Inner,2:5)= -1;
                                                iter_solv_Inner= iter_solv_Inner-1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        %
                        for i=indx_calc_delta_outer        % i_slct_delta_outer
                            if  termint_time_indc==1
                                break
                            end
                            %
                            for j=indx_calc_pi_outer
                                fultime= toc(teek0);
                                solution_time= fultime-Modeling_Time+solving_time;
                                if solution_time>= termint_time
                                    termint_time_indc=1;
                                    break
                                end

                                indcat_forbiden=0;
                                if size(unique([Pi(j).tree, Delta(i).forbiden_pi]),2) < size([Pi(j).tree, Delta(i).forbiden_pi],2)
                                    indcat_forbiden=1;
                                end
                                if indcat_forbiden == 0
                                    fultime= toc(teek0);
                                    solution_time= fultime-Modeling_Time+solving_time;
                                    if solution_time>= termint_time
                                        termint_time_indc=1;
                                        break
                                    else
                                        termint_time_func= termint_time-solution_time;
                                    end
                                    iter_solv_Outer= iter_solv_Outer+1;
                                    teek=tic;
                                    [Objective(fuliter).outer(iter_solv_Outer,1), deter_time, solving_func_time, sol_status_out, Objective(fuliter).outer_var(iter_solv_Outer)] = CLSC_func (Delta(i).outcornerU, Pi(j).outcornerU, Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
                                    Modeling_Time=Modeling_Time+toc(teek);
                                    Objective(fuliter).outer(iter_solv_Outer,2:3)=[i, j];
                                    solving_deter_time= solving_deter_time +deter_time;
                                    deter_time= 0;
                                    solving_time= solving_time + solving_func_time;
                                    solving_func_time= 0;
                                    if sol_status_out ~= 0
                                        disp(['sol_status_out= ', num2str(sol_status_out)]);
                                        disp(['Solving Time= ', num2str(solving_time)]);
                                        disp(['Solution Time= ', num2str(solution_time)]);
                                        disp(['Deterministic Time= ', num2str(solving_deter_time)]);
                                        fultime= toc(teek0);
                                        disp(['Full Time=  ', num2str(fultime)]);
                                        disp(['Modeling Time= ', num2str(Modeling_Time)]);
                                        Objective(fuliter).outer(iter_solv_Outer,1)= -Inf;
                                        Objective(fuliter).outer(iter_solv_Outer,2:3)= -1;
                                        iter_solv_Outer= iter_solv_Outer-1;
                                    end
                                end
                            end
                        end
                    end
                end
                %% End of Calculate Obj func
                %% Difference
                ind_maxout_delta=0;
                ind_maxout_pi=0;
                disp('----------------------------------');
                [sorted_in, sorted_row_in]= sort(Objective(fuliter).inner(:,1),'descend');
                maxin(fuliter)=sorted_in(1);
                rowmaxin=sorted_row_in(1);
                disp(['Maxin=  ',num2str(maxin(fuliter))]);
                disp(['Row Maxin=  ',num2str(rowmaxin)]);
                best_clmn_indc_delta_scenario=Objective(fuliter).inner(rowmaxin,4);
                best_clmn_indc_pi_scenario=Objective(fuliter).inner(rowmaxin,5);
                if maxin(fuliter) >= Best.Inner
                    Best.Inner= maxin(fuliter);
                    Best.Inner_iter_row=[fuliter, rowmaxin];
                    if best_clmn_indc_delta_scenario == 0
                        Best.incornerU_delta= Delta(Objective(fuliter).inner(rowmaxin,2)).incornerU;
                    else
                        Best.incornerU_delta= Delta(Objective(fuliter).inner(rowmaxin,2)).scenario(:,best_clmn_indc_delta_scenario);
                    end
                    if best_clmn_indc_pi_scenario == 0
                        Best.incornerU_pi= Pi(Objective(fuliter).inner(rowmaxin,3)).incornerU;
                    else
                        Best.incornerU_pi= Pi(Objective(fuliter).inner(rowmaxin,3)).scenario(:,best_clmn_indc_pi_scenario);
                    end
                    disp('new Best maxin');
                else
                    disp('NO new Best maxin');
                end
                %
                [maxout(fuliter), rowmaxout]= max(Objective(fuliter).outer(:,1));
                disp(['Maxout= ',num2str(maxout(fuliter))]);
                disp(['Row Maxout=  ',num2str(rowmaxout)]);
                %
                for counter_outer_obj= 1:iter_solv_Outer
                    if Objective(fuliter).outer(counter_outer_obj,1) < Best.Inner
                        Delta(Objective(fuliter).outer(counter_outer_obj,2)).forbiden_pi= [Delta(Objective(fuliter).outer(counter_outer_obj,2)).forbiden_pi, Objective(fuliter).outer(counter_outer_obj,3)];
                    end
                end
                %
                if maxout(fuliter) > Best.Inner
                    ind_maxout_delta= Objective(fuliter).outer(rowmaxout,2);
                    ind_maxout_pi= Objective(fuliter).outer(rowmaxout,3);
                    if maxout(fuliter) <= Best.Outer
                        Best.Outer= maxout(fuliter);
                        Best.Outer_iter_row=[fuliter, rowmaxout];
                        Best.outcornerU_delta= Delta(ind_maxout_delta).outcornerU;
                        Best.outcornerU_pi= Pi(ind_maxout_pi).outcornerU;
                        disp('new Best maxout');
                    else
                        if Best.Inner > Best.Outer % Here best inner has been improved and takes larger value than best outer (maxout > best inner).
                            Best.Outer= maxout(fuliter);
                            Best.Outer_iter_row=[fuliter, rowmaxout];
                            Best.outcornerU_delta= Delta(ind_maxout_delta).outcornerU;
                            Best.outcornerU_pi= Pi(ind_maxout_pi).outcornerU;
                            disp('new Best maxout (Best outer is elevated)');
                        else
                            disp('NO new Best maxout');
                        end
                    end
                    disp('----------------------------------');
                    % delat
                    if ismember(ind_maxout_delta, Delta(fuliter).chosen_calc(1, 1:siz_slct_delta(fuliter,1))) == 0  % yes >> maxout from final outers
                        disp(['Best Outer: Delta (', num2str(ind_maxout_delta),'): From final leftovers']);
                        ind_maxout_delta=0; % >> End of Delta wave
                    else
                        Delta(fuliter+1).chosen_split= ind_maxout_delta;
                        disp(['Best Outer: Delta (', num2str(ind_maxout_delta),') -> next iter']);
                    end
                    %pi
                    if ismember(ind_maxout_pi, Pi(fuliter).chosen_calc(1, 1:siz_slct_pi(fuliter,1))) == 0  % yes >> maxout from final outers
                        disp(['Best Outer: pi (', num2str(ind_maxout_pi),'): From final leftovers']);
                        ind_maxout_pi=0;  % >> End of Pi wave
                    else
                        Pi(fuliter+1).chosen_split= ind_maxout_pi;
                        disp(['Best Outer: pi (', num2str(ind_maxout_pi),') -> next iter']);
                    end
                else
                    if Best.Inner > Best.Outer % Here best inner has been improved and takes larger value than best outer (maxout < best inner).
                        Best.Outer= Inf;
                        Best.Outer_iter_row= 0;
                        Best.outcornerU_delta= 0;
                        Best.outcornerU_pi= 0;
                        disp('No new Best outer but Best outer is elevated to continue the process');
                    else
                        disp('maxout < maxin');    % >> End of Delta and Pi waves
                    end                    
                end
                disp('----------------------------------');
                %% getting reserve
                out_in_differnc(fuliter)= Best.Outer - Best.Inner;
                Best.Differenc= out_in_differnc(fuliter);
                if Best.Differenc < check_Best_diff
                    check_Best_diff= Best.Differenc;
                    ind_check_Best_diff=1;
                else
                    ind_check_Best_diff= ind_check_Best_diff +1;
                end
                if ind_check_Best_diff >= 2*D
                    disp(['"Min difference has not been changed for (',num2str(ind_check_Best_diff),') iteration >>> stop getting reserve"']);
                    disp('----------------------------------');
                else
                    disp('"Getting reserve"');
                    disp('----------------------------------');
                    % delta reserve
                    if ind_maxout_delta == 0
                        disp('"End of Delta wave, no reserve is gotten"');
                        disp('----------------------------------');
                    else
                        for counter_Objective_inner= 1:iter_solv_Inner
                            if Objective(fuliter).inner(sorted_row_in(counter_Objective_inner),4) == 0   % yes >> it is not a scenario
                                ind_maxin_delta= Objective(fuliter).inner(sorted_row_in(counter_Objective_inner),2);
                                if ind_maxin_delta ~= ind_maxout_delta
                                    if ismember(ind_maxin_delta, Delta(fuliter).chosen_calc(1, 1:siz_slct_delta(fuliter,1))) == 1  % yes >> maxin from this iteration
                                        reserv_delta(1, end+1)= ind_maxin_delta;
                                        break
                                    end
                                end
                            end
                        end
                    end
                    % pi reserve
                    if ind_maxout_pi == 0
                        disp('"End of Pi wave, no reserve is gotten"');
                        disp('----------------------------------');
                    else
                        for counter_Objective_inner= 1:iter_solv_Inner
                            if Objective(fuliter).inner(sorted_row_in(counter_Objective_inner),5) == 0   % yes >> it is not a scenario
                                ind_maxin_pi= Objective(fuliter).inner(sorted_row_in(counter_Objective_inner),3);
                                if ind_maxin_pi ~= ind_maxout_pi
                                    if ismember(ind_maxin_pi, Pi(fuliter).chosen_calc(1, 1:siz_slct_pi(fuliter,1))) == 1  % yes >> maxin from this iteration
                                        reserv_pi(1, end+1)= ind_maxin_pi;
                                        break
                                    end
                                end
                            end
                        end
                    end
                end
                %% wiping
                if maxin(fuliter) == Best.Inner
                    % delta wipe
                    ind_maxin_delta= Objective(fuliter).inner(sorted_row_in(1),2);
                    if ismember(ind_maxin_delta, Delta(fuliter).chosen_calc) == 1
                        if best_clmn_indc_delta_scenario == 0
                            disp(['Best Inner: Delta [',num2str(ind_maxin_delta),']']);
                            siz_slct_delta(fuliter,1)=1;
                            Delta(fuliter).chosen_calc=ind_maxin_delta;
                            siz_slct_delta_scenario(fuliter,1)=0;
                            Delta(fuliter).chosen_scenario_colmn=0;
                        else
                            disp(['Best Inner: Delta [',num2str(ind_maxin_delta),'] scenario(',num2str(best_clmn_indc_delta_scenario),')']);
                            siz_slct_delta(fuliter,1)=0;
                            Delta(fuliter).chosen_calc=ind_maxin_delta;
                            siz_slct_delta_scenario(fuliter,1)=1;
                            Delta(fuliter).chosen_scenario_colmn=best_clmn_indc_delta_scenario;
                        end
                    else
                        if best_clmn_indc_delta_scenario == 0
                            disp(['Best Inner: Delta [',num2str(ind_maxin_delta),']: From previuos']);
                        else
                            disp(['Best Inner: Delta [',num2str(ind_maxin_delta),'] scenario(',num2str(best_clmn_indc_delta_scenario),'): From previuos']);
                        end
                        siz_slct_delta(fuliter,1)=0;
                        Delta(fuliter).chosen_calc=double.empty(1,0);
                        siz_slct_delta_scenario(fuliter,1)=0;
                        Delta(fuliter).chosen_scenario_colmn=double.empty(1,0);
                    end
                    % pi wipe
                    ind_maxin_pi= Objective(fuliter).inner(sorted_row_in(1),3);
                    if ismember(ind_maxin_pi, Pi(fuliter).chosen_calc) == 1
                        if best_clmn_indc_pi_scenario == 0
                            disp(['Best Inner: Pi [',num2str(ind_maxin_pi),']']);
                            siz_slct_pi(fuliter,1)=1;
                            Pi(fuliter).chosen_calc=ind_maxin_pi;
                            siz_slct_pi_scenario(fuliter,1)=0;
                            Pi(fuliter).chosen_scenario_colmn=0;
                        else
                            disp(['Best Inner: Pi [',num2str(ind_maxin_pi),'] scenario(',num2str(best_clmn_indc_pi_scenario),')']);
                            siz_slct_pi(fuliter,1)=0;
                            Pi(fuliter).chosen_calc=ind_maxin_pi;
                            siz_slct_pi_scenario(fuliter,1)=1;
                            Pi(fuliter).chosen_scenario_colmn=best_clmn_indc_pi_scenario;
                        end
                    else
                        if best_clmn_indc_pi_scenario == 0
                            disp(['Best Inner: Pi [',num2str(ind_maxin_pi),']: From previuos']);
                        else
                            disp(['Best Inner: Pi [',num2str(ind_maxin_pi),'] scenario(',num2str(best_clmn_indc_pi_scenario),'): From previuos']);
                        end
                        siz_slct_pi(fuliter,1)=0;
                        Pi(fuliter).chosen_calc=double.empty(1,0);
                        siz_slct_pi_scenario(fuliter,1)=0;
                        Pi(fuliter).chosen_scenario_colmn=double.empty(1,0);
                    end
                else
                    siz_slct_delta(fuliter,1)=0;
                    Delta(fuliter).chosen_calc=double.empty(1,0);
                    siz_slct_delta_scenario(fuliter,1)=0;
                    Delta(fuliter).chosen_scenario_colmn=double.empty(1,0);
                    %
                    siz_slct_pi(fuliter,1)=0;
                    Pi(fuliter).chosen_calc=double.empty(1,0);
                    siz_slct_pi_scenario(fuliter,1)=0;
                    Pi(fuliter).chosen_scenario_colmn=double.empty(1,0);
                    disp('Best Inner: nothing better was found ');
                end
                disp('----------------------------------');
                disp(['Reserve delta:[', num2str(reserv_delta) ,']']);
                disp(['Reserve pi:   [', num2str(reserv_pi) ,']']);
                if size(Delta(fuliter+1).chosen_split,2) ==0
                    if size(reserv_delta,2) ~= 0
                        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
                        Delta(fuliter+1).chosen_split= reserv_delta(1); % FIFO
                        disp(['chosen split delta changed to Reserve: [', num2str(reserv_delta(1)) ,']']);
                        reserv_delta(1)=[];
                    else
                        siz_old_acceptbl_delta=0;
                    end
                end
                if size(Pi(fuliter+1).chosen_split,2) ==0
                    if size(reserv_pi,2) ~= 0
                        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
                        Pi(fuliter+1).chosen_split= reserv_pi(1); % FIFO
                        disp(['chosen split pi changed to Reserve: [', num2str(reserv_pi(1)) ,']']);
                        reserv_pi(1)=[];
                    else
                        siz_old_acceptbl_pi=0;
                    end
                end
                fuliter_solv_Inner= fuliter_solv_Inner +iter_solv_Inner;
                fuliter_solv_Outer= fuliter_solv_Outer +iter_solv_Outer;
                %
                if out_in_differnc(fuliter) <= epsilon
                    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
                    disp('"Termination: Maxout and Maxin differ reached"');
                    fultime= toc(teek0);
                    break
                end
            end
            disp('----------------------------------');
            disp(['Total maxin=  ',num2str(Best.Inner)]);
            disp(['Total maxout= ',num2str(Best.Outer)]);
            disp(['Min difference= ', num2str(Best.Differenc)]);
            disp('----------------------------------');
            if termint_time_indc==1
                disp('Termination: Time termination')
                break
            end
            %% End of Difference and calculation
            fultime= toc(teek0);
            solution_time= fultime-Modeling_Time+solving_time;
            disp(['Full Time=  ', num2str(fultime)]);
            disp(['Modeling Time= ', num2str(Modeling_Time)]);
            disp(['Solution Time= ', num2str(solution_time)]);
            disp(['Deterministic Time= ', num2str(solving_deter_time)]);
            disp(['Alg24 solver Time= ', num2str(solving_Alg24_time)]);
            disp(['Alg7 time: ',num2str(Alg7_time)]);
            disp(['dif_solVSalg24_time: ',num2str(dif_solVSalg24_time)]);
            disp(['deter+alg24+dif_solVSalg24_time: ',num2str(solving_deter_time+solving_Alg24_time+dif_solVSalg24_time)]);
            disp(['Number of inner solving: ', num2str(fuliter_solv_Inner)])
            disp(['Number of outer solving: ', num2str(fuliter_solv_Outer)])
            disp(['Number of All solving: ', num2str(fuliter_solv_Inner+fuliter_solv_Outer)])
            %
            if siz_old_acceptbl_delta ==0 && siz_old_acceptbl_pi ==0
                disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
                disp('All Deltas and Pis have been checked');
                break
            elseif siz_old_acceptbl_delta ==0
                disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
                disp('All Deltas have been checked');
            elseif siz_old_acceptbl_pi ==0
                disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
                disp('All Pis have been checked');
            end
        end
        %% End of Main loop
        %% Results
        disp('###################################');
        if  termint_time_indc==2
            disp(['Final completed Iteration (', num2str(fuliter-1),') Results']);
        else
            disp(['Final completed Iteration (', num2str(fuliter),') Results']);
        end
        disp('***********************************');
        solution_time= fultime-Modeling_Time+solving_time;  % ``solving_deter_time" is a part of ``solving_time" is a part of ``Modeling_Time" is a part of ``fultime"
        disp(['Full Time=  ', num2str(fultime)]);
        disp(['Modeling Time= ', num2str(Modeling_Time)]);
        disp(['Solution Time= ', num2str(solution_time)]);
        disp(['Deterministic Time= ', num2str(solving_deter_time)]);
        disp(['Alg24 solver Time= ', num2str(solving_Alg24_time)]);
        disp(['Alg7 time: ',num2str(Alg7_time)]);
        disp(['dif_solVSalg24_time: ',num2str(dif_solVSalg24_time)]);
        disp(['deter+alg24+dif_solVSalg24_time: ',num2str(solving_deter_time+solving_Alg24_time+dif_solVSalg24_time)]);
        disp(' ');
        disp(['All differs: [', num2str(out_in_differnc) ,']'])
        disp(' ');
        disp(['Best_Inner= ', num2str(Best.Inner)]);
        disp(['Best_outer= ', num2str(Best.Outer)]);
        disp(['Min difference= ', num2str(Best.Differenc)]);
        Best.error= (Best.Differenc/Best.Inner)*100;
        disp(['Error= ', num2str(Best.error), '%']);
        disp(' ');
        disp(['Best_incornerU_delta: [', num2str(Best.incornerU_delta') ,']'])
        disp(['Best_incornerU_pi: [', num2str(Best.incornerU_pi') ,']'])
        disp(['Best_outcornerU_delta: [', num2str(Best.outcornerU_delta') ,']'])
        disp(['Best_outcornerU_pi: [', num2str(Best.outcornerU_pi') ,']'])
        disp(' ');
        disp(['Number of inner solving: ', num2str(fuliter_solv_Inner)])
        disp(['Number of outer solving: ', num2str(fuliter_solv_Outer)])
        disp(['Number of All solving: ', num2str(fuliter_solv_Inner+fuliter_solv_Outer)])
        disp('*********************************');
        %% End of Results
        diary off
        save([pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Heuristic_v12\',MatName,'_H_v12_',num2str(detergap),'.mat'],'Best','Delta','Pi','Objective','maxin','maxout','PoolscenarioU_delta','PoolscenarioU_pi','sol_status_in','sol_status_out','siz_scenarioU_delta','siz_scenarioU_pi','siz_slct_delta','siz_slct_delta_scenario','siz_slct_pi','siz_slct_pi_scenario','siz_Total_delta','siz_Total_pi','fultime','Modeling_Time','solving_time','solution_time','solving_deter_time','solving_Alg24_time','termint_time_func','termint_time_indc');
        zip([pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Heuristic_v12\',MatName,'_H_v12_',num2str(detergap)],{[pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Heuristic_v12\', MatName,'_H_v12_',num2str(detergap),'.mat'],[pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Heuristic_v12\',MatName, '_H_v12_', num2str(detergap),'.txt']});
        delete([pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Heuristic_v12\',MatName,'_H_v12_',num2str(detergap),'.mat']);
        delete([pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Heuristic_v12\',MatName,'_H_v12_',num2str(detergap),'.txt']);
        close all;
        clearvars -except run_iter
        clear functions
        clear global
        clc;
    end
end
%
% Initial inventory level of each brewery shuold be lower than total filled bottles demand