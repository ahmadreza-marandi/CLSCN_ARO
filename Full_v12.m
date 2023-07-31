%% Main_Alg
close all;
clear all;
clc;
format short
diary off
warning('off')
%
for run_iter= 7
    for ins_iter= 2:5
        MatName= ['ins_', num2str(run_iter), '_iter_', num2str(ins_iter), '_V8'];
        if run_iter==4 || run_iter==7
            detergap= .1;
        elseif run_iter==10
            detergap= .15;
        else
            detergap= 10^(-4);
        end
        DiaryName= [MatName, '_Full_v12_', num2str(detergap),'.txt'];
        diary([pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Full_v12\', DiaryName])
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
        termint_time=400*D;    % On solution time
        dif_incorner_U_L=2;   % decimal acceptance : 2   % delta and pi spliting acceptance : 0.01
        epsilon=.01;
        diffVol_delta=.001;
        diffVol_pi=.001;
        disp(['Uper-Lower Termination=  ', num2str(epsilon)]);
        disp(['Delta Vol out-in Ter.=   ', num2str(diffVol_delta)]);
        disp(['Pi Vol out-in Ter.=      ', num2str(diffVol_pi)]);
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
        dif_out_in_delta(1,1)= Delta(1).Volout - Delta(1).Volin;
        % first delta selection
        siz_slct_delta(fuliter,1)=1;                                                      % #size delta select (solve)
        siz_scenarioU_delta(1,1)= size(Delta(1).scenario,2);
        siz_slct_delta_scenario(fuliter,1)=  siz_scenarioU_delta(1,1);               % #size delta scenario select (solve)
        Delta(fuliter).chosen_calc= ones(1,siz_slct_delta_scenario(fuliter,1)+1);  % One delta wil be selected
        Delta(fuliter).chosen_scenario_colmn= [0, 1:siz_slct_delta_scenario(fuliter,1)];  % column(s) of delta scenario
        Delta(fuliter+1).chosen_split= 1;
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
        dif_out_in_pi(1,1)= Pi(1).Volout- Pi(1).Volin;
        % first pi selection
        siz_slct_pi(fuliter,1)=1;                                                        % #size pi select (solve)
        siz_scenarioU_pi(1,1)= size(Pi(1).scenario,2);
        siz_slct_pi_scenario(fuliter,1)= siz_scenarioU_pi(1,1);                    % #size pi scenario select (solve)
        Pi(fuliter).chosen_calc= ones(1,siz_slct_pi_scenario(fuliter,1)+1);  % One pi wil be selected
        Pi(fuliter).chosen_scenario_colmn= [0, 1:siz_slct_pi_scenario(fuliter,1)];           % column(s) of pi scenario
        Pi(fuliter+1).chosen_split= 1;
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
                        [Objective(fuliter).inner(iter_solv_Inner,1), deter_time, solving_func_time, sol_status_in, Objective(fuliter).inner_var(iter_solv_Inner)] = CLSC_func (Delta(1).incornerU, Pi(1).scenario(:,j-1), Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
                        Objective(fuliter).inner(iter_solv_Inner,2:5)=[1, 1, 0, j-1];    % [ delta, pi, scenario_delta, scenario_pi ]

                    end
                else
                    if j == 1 % delta scenarios, pi inners
                        [Objective(fuliter).inner(iter_solv_Inner,1), deter_time, solving_func_time, sol_status_in, Objective(fuliter).inner_var(iter_solv_Inner)] = CLSC_func (Delta(1).scenario(:,i-1), Pi(1).incornerU, Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
                        Objective(fuliter).inner(iter_solv_Inner,2:5)=[1, 1, i-1, 0];
                    else % both scenarios
                        [Objective(fuliter).inner(iter_solv_Inner,1), deter_time, solving_func_time, sol_status_in, Objective(fuliter).inner_var(iter_solv_Inner)] = CLSC_func (Delta(1).scenario(:,i-1), Pi(1).scenario(:,j-1), Obj_func, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap); % Alg CLSC
                        Objective(fuliter).inner(iter_solv_Inner,2:5)=[1, 1, i-1, j-1];    % [ delta, pi, scenario_delta, scenario_pi ]
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
                disp(['sol_status_out= ', num2str(sol_status_out)]);
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
        elseif best_clmn_indc_delta_scenario ~= -1
            disp(['Best Inner: delta [1] scenario(',num2str(best_clmn_indc_delta_scenario),')'])
            Best.incornerU_delta= Delta(1).scenario(:,best_clmn_indc_delta_scenario);
        end
        if best_clmn_indc_pi_scenario == 0
            disp('Best Inner: pi [1]')
            Best.incornerU_pi= Pi(1).incornerU;
        elseif best_clmn_indc_pi_scenario ~= -1
            disp(['Best Inner: pi [1] scenario(',num2str(best_clmn_indc_pi_scenario),')'])
            Best.incornerU_pi= Pi(1).scenario(:,best_clmn_indc_pi_scenario);
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
        disp(['Number of inner solving: ', num2str(fuliter_solv_Inner)])
        disp(['Number of outer solving: ', num2str(fuliter_solv_Outer)])
        disp(['Number of All solving: ', num2str(fuliter_solv_Inner+fuliter_solv_Outer)])
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
            Delta(fuliter).chosen_scenario_colmn= double.empty(1,0);
            Pi(fuliter).chosen_scenario_colmn= double.empty(1,0);
            siz_acceptbl_delta= 0;
            siz_acceptbl_pi= 0;
            siz_slct_delta(fuliter,1)= 0;
            siz_slct_pi(fuliter,1)= 0;
            siz_slct_delta_scenario(fuliter,1)=0;
            siz_slct_pi_scenario(fuliter,1)=0;
            %% Delta iteration
            disp('==============(Delta)==============');
            if siz_old_acceptbl_delta~=0
                disp(['siz_Total_delta until now= {', num2str(siz_Total_delta),'}']);
                disp('- - - - - - - - - - -');
                for iter_split_delta= Delta(fuliter).chosen_split
                    if  termint_time_indc==2
                        break
                    end
                    siznewA_delta=0;
                    if dif_out_in_delta(iter_split_delta,1) > diffVol_delta
                        teek=tic;
                        [newA, newb, solving_func_time]= Alg7_func (Delta(iter_split_delta).A, Delta(iter_split_delta).b, D, Delta(iter_split_delta).incornerL, Delta(iter_split_delta).incornerU, prelambda, lambda); % Alg 7
                        Modeling_Time=Modeling_Time+toc(teek);
                        solving_time= solving_time + solving_func_time;
                        solving_func_time= 0;
                        siznewA_delta= size(newA,2);
                        for i= 1:siznewA_delta
                            fultime= toc(teek0);
                            Modeling_Time=Modeling_Time;
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
                            solving_alg24_func_time= 0;
                            solving_time= solving_time + solving_func_time;
                            solving_func_time= 0;
                            if ismember(0, round(Delta(siz_Total_delta).incornerU - Delta(siz_Total_delta).incornerL,dif_incorner_U_L)) == 0
                                siz_acceptbl_delta= siz_acceptbl_delta+ 1;
                                siz_scenarioU_delta(siz_Total_delta,1)=size( Delta(siz_Total_delta).scenario,2);
                                Delta(siz_Total_delta).forbiden_pi=Delta(iter_split_delta).forbiden_pi;
                                %
                                Delta(siz_Total_delta).A= newA(i).part;
                                Delta(siz_Total_delta).b= newb(i).part;
                                dif_out_in_delta(siz_Total_delta,1)= Delta(siz_Total_delta).Volout-Delta(siz_Total_delta).Volin; % it's on Totalsize delta
                            else
                                siz_Total_delta= siz_Total_delta- 1;
                                disp('"One new Delat polytope rejected"');
                                siznewA_delta=siznewA_delta- 1;
                            end
                        end
                        disp(['Delta(', num2str(iter_split_delta),'): "', num2str(dif_out_in_delta(iter_split_delta,1)),'", +',num2str(siznewA_delta),' = ', num2str(siz_Total_delta)]);
                        if siznewA_delta == 0
                            final_Outer_delta= final_Outer_delta+1;
                            indx_final_outer_delta(1,final_Outer_delta)= iter_split_delta;  % for OuterU Delta belongs to bad shaped polyttops
                            disp(['Delta(', num2str(iter_split_delta),'): "', num2str(dif_out_in_delta(iter_split_delta,1)),'" -> final_Outer_delta (',num2str(final_Outer_delta),') [Because of shape]']);
                        end
                    else
                        final_Outer_delta= final_Outer_delta+1;
                        indx_final_outer_delta(1,final_Outer_delta)= iter_split_delta;  % for OuterU Delta belongs to small polyttops
                        disp(['Delta(', num2str(iter_split_delta),'): "', num2str(dif_out_in_delta(iter_split_delta,1)),'" -> final_Outer_delta (',num2str(final_Outer_delta),') [Because of size]']);
                    end
                end
                if  termint_time_indc==2
                    break
                end
                siz_old_acceptbl_delta= siz_acceptbl_delta; % updating sizeoldA_delta for next Step
                if siz_old_acceptbl_delta ~= 0
                    siz_slct_delta(fuliter,1)= siz_old_acceptbl_delta;        % #size delta select (solve)
                    Delta(fuliter).chosen_calc= siz_Total_delta-siz_old_acceptbl_delta+1:siz_Total_delta;
                    Delta(fuliter).chosen_scenario_colmn= zeros(1, siz_slct_delta(fuliter,1));
                    disp('- - - - - - - - - - -');
                    % delta_scenario_choose
                    siz_slct_delta_scenario(fuliter,1)= sum(siz_scenarioU_delta(Delta(fuliter).chosen_calc));         % #size delta scenario select (solve)
                    for i=1:siz_slct_delta(fuliter,1)
                        Delta(fuliter).chosen_calc(1,end+1:end+siz_scenarioU_delta(Delta(fuliter).chosen_calc(i)))= Delta(fuliter).chosen_calc(i);
                        Delta(fuliter).chosen_scenario_colmn(1,end+1:end+siz_scenarioU_delta(Delta(fuliter).chosen_calc(i)))= (1:siz_scenarioU_delta(Delta(fuliter).chosen_calc(i)));
                    end
                    disp(['chosen delta to calculat: [', num2str(Delta(fuliter).chosen_calc) ,']'])
                    disp(['scenario: [', num2str(Delta(fuliter).chosen_scenario_colmn) ,']'])
                else
                    disp('* No new Delta from split');
                end
            else
                disp('* No new Delta from split');
            end
            %% End of Delta iteration
            %% pi iteration
            disp('==============(Pi)==============');
            if siz_old_acceptbl_pi~=0
                disp(['siz_Total_pi until now= {', num2str(siz_Total_pi),'}']);
                disp('- - - - - - - - - - -');
                for iter_split_pi= Pi(fuliter).chosen_split
                    if  termint_time_indc==2
                        break
                    end
                    siznewA_pi=0;
                    if dif_out_in_pi(iter_split_pi,1) > diffVol_pi
                        teek=tic;
                        [newA, newb, solving_func_time]= Alg7_func (Pi(iter_split_pi).A, Pi(iter_split_pi).b, D, Pi(iter_split_pi).incornerL, Pi(iter_split_pi).incornerU, prelambda, lambda); % Alg 7
                        Modeling_Time=Modeling_Time+toc(teek);
                        solving_time= solving_time + solving_func_time;
                        solving_func_time= 0;
                        siznewA_pi= size(newA,2);
                        for i= 1:siznewA_pi
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
                            solving_alg24_func_time= 0;
                            solving_time= solving_time + solving_func_time;
                            solving_func_time= 0;
                            if ismember(0, round(Pi(siz_Total_pi).incornerU - Pi(siz_Total_pi).incornerL,dif_incorner_U_L)) == 0
                                siz_acceptbl_pi= siz_acceptbl_pi+ 1;
                                siz_scenarioU_pi(siz_Total_pi,1)=size( Pi(siz_Total_pi).scenario,2);
                                Pi(siz_Total_pi).tree=[Pi(iter_split_pi).tree, siz_Total_pi];

                                Pi(siz_Total_pi).A= newA(i).part;
                                Pi(siz_Total_pi).b= newb(i).part;
                                dif_out_in_pi(siz_Total_pi,1)= Pi(siz_Total_pi).Volout-Pi(siz_Total_pi).Volin; % it's on Totalsize pi
                            else
                                siz_Total_pi= siz_Total_pi- 1;
                                disp('"One new Delat polytope rejected"');
                                siznewA_pi=siznewA_pi- 1;
                            end
                        end
                        disp(['Pi(', num2str(iter_split_pi),'): "', num2str(dif_out_in_pi(iter_split_pi,1)),'", +',num2str(siznewA_pi),' = ', num2str(siz_Total_pi)]);
                        if siznewA_pi == 0
                            final_Outer_pi= final_Outer_pi+1;
                            indx_final_outer_pi(1,final_Outer_pi)= iter_split_pi;  % for OuterU Pi belongs to bad shaped polyttops
                            disp(['Pi(', num2str(iter_split_pi),'): "', num2str(dif_out_in_pi(iter_split_pi,1)),'" -> final_Outer_pi (',num2str(final_Outer_pi),') [Because of shape]']);
                        end
                    else
                        final_Outer_pi= final_Outer_pi+1;
                        indx_final_outer_pi(1,final_Outer_pi)= iter_split_pi;  % for OuterU Pi belongs to small polyttops
                        disp(['Pi(', num2str(iter_split_pi),'): "', num2str(dif_out_in_pi(iter_split_pi,1)),'" -> final_Outer_pi (',num2str(final_Outer_pi),') [Because of size]']);
                    end
                end
                if  termint_time_indc==2
                    break
                end
                siz_old_acceptbl_pi= siz_acceptbl_pi; % updating sizeoldA_pi for next Step
                if siz_old_acceptbl_pi ~= 0
                    siz_slct_pi(fuliter,1)= siz_old_acceptbl_pi;        % #size pi select (solve)
                    Pi(fuliter).chosen_calc= siz_Total_pi-siz_old_acceptbl_pi+1:siz_Total_pi;
                    Pi(fuliter).chosen_scenario_colmn= zeros(1, siz_slct_pi(fuliter,1));
                    disp('- - - - - - - - - - -');
                    % pi_scenario_choose
                    siz_slct_pi_scenario(fuliter,1)= sum(siz_scenarioU_pi(Pi(fuliter).chosen_calc));         % #size pi scenario select (solve)
                    for i=1:siz_slct_pi(fuliter,1)
                        Pi(fuliter).chosen_calc(1,end+1:end+siz_scenarioU_pi(Pi(fuliter).chosen_calc(i)))= Pi(fuliter).chosen_calc(i);
                        Pi(fuliter).chosen_scenario_colmn(1,end+1:end+siz_scenarioU_pi(Pi(fuliter).chosen_calc(i)))= (1:siz_scenarioU_pi(Pi(fuliter).chosen_calc(i)));
                    end
                    disp(['chosen pi to calculat: [', num2str(Pi(fuliter).chosen_calc) ,']'])
                    disp(['scenario: [', num2str(Pi(fuliter).chosen_scenario_colmn) ,']'])
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
                disp('----------------------------------');
                [maxin(fuliter), rowmaxin]= max(Objective(fuliter).inner(:,1));
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
                if maxout(fuliter) ~= -Inf
                    ind_maxout_delta= Objective(fuliter).outer(rowmaxout,2);
                    ind_maxout_pi= Objective(fuliter).outer(rowmaxout,3);
                    if maxout(fuliter) <= Best.Outer
                        Best.Outer= maxout(fuliter);
                        Best.Outer_iter_row=[fuliter, rowmaxout];
                        Best.outcornerU_delta= Delta(ind_maxout_delta).outcornerU;
                        Best.outcornerU_pi= Pi(ind_maxout_pi).outcornerU;
                        disp('new Best maxout');
                    else
                        disp('NO new Best maxout');
                    end
                    disp('----------------------------------');
                    for counter_outer_obj= 1:iter_solv_Outer
                        if Objective(fuliter).outer(counter_outer_obj,1) < Best.Inner
                            Delta(Objective(fuliter).outer(counter_outer_obj,2)).forbiden_pi= [Delta(Objective(fuliter).outer(counter_outer_obj,2)).forbiden_pi, Objective(fuliter).outer(counter_outer_obj,3)];
                        end
                    end
                    % delat
                    if ismember(ind_maxout_delta, Delta(fuliter).chosen_calc(1, 1:siz_slct_delta(fuliter,1))) == 0 % yes >> maxout from final outers
                        disp(['Best Outer: Delta (', num2str(ind_maxout_delta),'): From final leftovers']);
                    else
                        disp(['Best Outer: Delta (', num2str(ind_maxout_delta),')']);
                    end
                    %pi
                    if ismember(ind_maxout_pi, Pi(fuliter).chosen_calc(1, 1:siz_slct_pi(fuliter,1))) == 0  % yes >> maxout from final outers
                        disp(['Best Outer: pi (', num2str(ind_maxout_pi),'): From final leftovers']);
                    else
                        disp(['Best Outer: pi (', num2str(ind_maxout_pi),')']);
                    end
                else
                    disp('Best Outer not found');
                end
                disp('----------------------------------');
                % delta inner
                ind_maxin_delta= Objective(fuliter).inner(rowmaxin,2);
                if ismember(ind_maxin_delta, Delta(fuliter).chosen_calc) == 1  % yes >> maxin from this iteration
                    if best_clmn_indc_delta_scenario == 0
                        disp(['Best Inner: Delta [',num2str(ind_maxin_delta),']']);
                    else
                        disp(['Best Inner: Delta [',num2str(ind_maxin_delta),'] scenario(',num2str(best_clmn_indc_delta_scenario),')']);
                    end
                else
                    if best_clmn_indc_delta_scenario == 0
                        disp(['Best Inner: Delta [',num2str(ind_maxin_delta),']: From previuos']);
                    elseif best_clmn_indc_delta_scenario ~= -1
                        disp(['Best Inner: Delta [',num2str(ind_maxin_delta),'] scenario(',num2str(best_clmn_indc_delta_scenario),'): From previuos']);
                    else
                        disp('Best Inner: Delta not found');
                    end
                end
                % pi inner
                ind_maxin_pi= Objective(fuliter).inner(rowmaxin,3);
                if ismember(ind_maxin_pi, Pi(fuliter).chosen_calc) == 1  % yes >> maxin from this iteration
                    if best_clmn_indc_pi_scenario == 0
                        disp(['Best Inner: Pi [',num2str(ind_maxin_pi),']']);
                    else
                        disp(['Best Inner: Pi [',num2str(ind_maxin_pi),'] scenario(',num2str(best_clmn_indc_pi_scenario),')']);
                    end
                else
                    if best_clmn_indc_pi_scenario == 0
                        disp(['Best Inner: Pi [',num2str(ind_maxin_pi),']: From previuos']);
                    elseif best_clmn_indc_pi_scenario ~= -1
                        disp(['Best Inner: Pi [',num2str(ind_maxin_pi),'] scenario(',num2str(best_clmn_indc_pi_scenario),'): From previuos']);
                    else
                        disp('Best Inner: Pi not found');
                    end
                end
                fuliter_solv_Inner= fuliter_solv_Inner +iter_solv_Inner;
                fuliter_solv_Outer= fuliter_solv_Outer +iter_solv_Outer;
                %
                out_in_differnc(fuliter)= Best.Outer - Best.Inner;
                Best.Differenc= out_in_differnc(fuliter);
                if out_in_differnc(fuliter) <= epsilon
                    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
                    disp('"Termination: Maxout and Maxin differ reached"');
                    fultime=toc(teek0);
                    break
                end
                Delta(fuliter+1).chosen_split= Delta(fuliter).chosen_calc(1:siz_slct_delta(fuliter,1));
                Pi(fuliter+1).chosen_split= Pi(fuliter).chosen_calc(1:siz_slct_pi(fuliter,1));
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
            disp(['Number of inner solving: ', num2str(fuliter_solv_Inner)])
            disp(['Number of outer solving: ', num2str(fuliter_solv_Outer)])
            disp(['Number of All solving: ', num2str(fuliter_solv_Inner+fuliter_solv_Outer)])
            %
            if (siz_old_acceptbl_delta ==0 && siz_old_acceptbl_pi ==0) || maxout(fuliter) == -Inf
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
        save([pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Full_v12\',MatName,'_Full_v12_',num2str(detergap),'.mat'],'Best','Delta','Pi','Objective','maxin','maxout','PoolscenarioU_delta','PoolscenarioU_pi','sol_status_in','sol_status_out','siz_scenarioU_delta','siz_scenarioU_pi','siz_slct_delta','siz_slct_delta_scenario','siz_slct_pi','siz_slct_pi_scenario','siz_Total_delta','siz_Total_pi','fultime','Modeling_Time','solving_time','solution_time','solving_deter_time','solving_Alg24_time','termint_time_func','termint_time_indc');
        zip([pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Full_v12\',MatName,'_Full_v12_',num2str(detergap)],{[pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Full_v12\', MatName,'_Full_v12_',num2str(detergap),'.mat'],[pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Full_v12\',MatName, '_Full_v12_', num2str(detergap),'.txt']});
        delete([pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Full_v12\',MatName,'_Full_v12_',num2str(detergap),'.mat']);
        delete([pwd, '\Runs_group\', '\ins_', num2str(run_iter), '\Full_v12\',MatName,'_Full_v12_',num2str(detergap),'.txt']);
        close all;
        clearvars -except run_iter
        clear functions
        clear global
        clc;
    end
end
%
% Initial inventory level of each brewery shuold be lower than total filled bottles demand