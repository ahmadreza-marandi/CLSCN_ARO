function [ outcornerU, incornerL, incornerU, scenario, scenario_0, solving_alg24_func_time, solving_func_time, Volin, Volout] = Alg24_func (A_ini ,b_ini, x, D, scenario_0)
%Inner and Outer calc with scenarios
%% Initial Parameters
b_ini= round(b_ini, 10);
size_scenario_0=size(scenario_0,2);
epsilon_func=0.0001;
outcornerL(1:D,1)=0;
outcornerU(1:D,1)=0;
scenario_prim= double.empty(D,0);
modeling_time=0;
solving_alg24_func_time=0;
teek0=tic;
%% Outer Box (Algorithm 4)
for i=1:D    % Minimization    
    teek=tic;
    obj= x(i);
    sol1=optimize(A_ini*x <= b_ini -epsilon_func, obj, sdpsettings('solver','cplex', 'verbose',0, 'cachesolvers',1));
    modeling_time=modeling_time+toc(teek);
    solving_alg24_func_time= solving_alg24_func_time+sol1.solvertime;   % solving deter time
    outcornerL(i,1)= value(x(i));
end
for i=1:D    % Maximization    
    teek=tic;
    obj= -x(i);
    sol2=optimize(A_ini*x <= b_ini , obj, sdpsettings('solver','cplex', 'verbose',0, 'cachesolvers',1));
    modeling_time=modeling_time+toc(teek);
    solving_alg24_func_time= solving_alg24_func_time+sol2.solvertime;
    outcornerU(i,1)= value(x(i));    
    scenario_prim= [scenario_prim, value(x)];
end
Volout= prod(outcornerU - outcornerL);
%% Inner Box
teek=tic;
Innerlambda = sdpvar(1);
Aplus= max(A_ini,0);
r= sdpvar(D,1);
innerP3= [A_ini*x + Aplus*r*Innerlambda <= b_ini; r>=0 ; Innerlambda>=0];
assign(r, ones(D,1));
assign(x, outcornerL );
assign(Innerlambda, min(mean([outcornerU,outcornerL],2) -outcornerL));
Objective3= - sum(r*Innerlambda);
sol3=optimize(innerP3, Objective3, sdpsettings('solver','FMINCON','verbose',0,'usex0',1, 'cachesolvers',1));
modeling_time=modeling_time+toc(teek);
solving_alg24_func_time= solving_alg24_func_time+sol3.solvertime;
incornerL=value(x);
incornerU=value(x) + value(r)*value(Innerlambda);
Volin= prod(incornerU - incornerL);
%%
scenario_prim= [scenario_prim, incornerU];
scenario_0_prim= [scenario_0, scenario_prim];
[~,Iv]=unique(round(scenario_0_prim',2),'rows');
if size(Iv,1) == size_scenario_0
    scenario= [];
else
    Iv=sort(Iv);
    Iv(1:size_scenario_0)=[];
    Iv=Iv-size_scenario_0;
    scenario_0=[scenario_0, scenario_prim(:,Iv)];
    Iv(end)=[];
    scenario= scenario_prim(:,Iv);
end
%%
tFinl=toc(teek0);
solving_func_time= tFinl-modeling_time+solving_alg24_func_time;
end