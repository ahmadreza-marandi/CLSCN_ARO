function [Obj, deter_time, solving_func_time, sol_status, val] = CLSC_func (cornerU_delta, cornerU_pi, Objective, Constraints, M, N, T, D, V, Iemp, Iful, f, e, purch, m, r, v, u, Q, X, Le, termint_time_func, detergap)
%% Parameters
delta=zeros(T,N);
delta(1:D)=cornerU_delta(1:D);
delta=[zeros(M,T);delta'];
pi=zeros(T,N);
pi(1:D)=cornerU_pi(1:D);
pi=[zeros(M,T);pi'];
%% Additional constraints
maxQ= max(Q(1,1,:,1));
for period=1:T
    Constraints = [Constraints;
        Iemp(M+1:V,period) - sum(e(:,:,period),2) + pi(M+1:V ,period)   <= Iemp(M+1:V,period+1);%Cons 4
        sum(f(:,:,period),2) >= delta(M+1:V ,period);%Cons 5
%         sum(sum(sum(X(M+1:V, 1:V, :, 1:period),2), 3), 4) >= ceil( sum(delta(M+1:V,1:period),2) ./ maxQ );  %25_paper (remove: because in our study, a retailer can be visited just by one vehicle, but in the base study (Qiu 2018) there is no limitation)
        sum(X(1:M, 1:V, :, 1:period),'all') >= ceil(sum(delta(M+1:V,1:period),[1,2]) ./ maxQ ) ];  %26_paper
%         sum(sum(sum(X(M+1:V, 1:V, :, 1:period),2), 3), 4) >= ceil((Iemp(M+1:V,1) + sum(pi(M+1:V,1:period),2) - Le(M+1:V) ) ./ min(maxQ,Le(M+1:V)) ) ];  %27_paper (remove: because in our study, a retailer can be visited just by one vehicle, but in the base study (Qiu 2018) there is no limitation)
        for s=0:period-1
            Constraints = [Constraints;
                sum(delta(M+1:V, period-s:period),2) .* (1- sum(sum(sum(X(M+1:V, :, :, 1:period), 2), 3), 4)) <= 0;  %33_paper
                Le(M+1:V) - sum(pi(M+1:V, period-s:period),2) .* (1- sum(sum(sum(X(M+1:V, :, :, 1:period), 2), 3), 4)) >= Iemp(M+1:V,period-s) ];  %34_paper
        end
end
%% Solution
sol= optimize(Constraints, Objective, sdpsettings('solver','cplex','verbose',0,'cplex.mip.tolerances.mipgap', detergap, 'cplex.timelimit', termint_time_func, 'cachesolvers',1));
deter_time= sol.solvertime;
Obj= value(Objective);
solving_func_time= deter_time;
sol_status= sol.problem;

val.purch= value(purch);
val.m= value(m);
val.r= value(r);
val.Iful= value(Iful);
val.Iemp= value(Iemp);
val.f= value(f);
val.e= value(e);
val.v= value(v);
val.u= value(u);
val.X= value(X);
val.delta=delta;
val.pi=pi;
end