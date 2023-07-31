function [Obj, deter_time, solving_func_time, sol_status, val] = RSOME_func (A_D, b_D, A_P, b_P, M, N, K, T, V, base_Q, base_FixC, base_C, FilCap, Cpur, Chyg, Cins, Cdis, He, Hf, Lf, Le, alpha, Iful0, Iemp0, termint_time_func, detergap)
b_D= round(b_D, 10);
b_P= round(b_P, 10);
model=rsome('RO_CLSC_Total');
model.Param.display=0;
model.Param.timelimit= termint_time_func;
model.Param.mipgap = detergap;
%
zeta= model.random(2*T*N);
delta=zeta(1:T*N);
pi=zeta(T*N+1:2*T*N);
%% Converting
C= reshape(permute(base_C,[1 3 2 4]),K*V,T*V);
FixC= reshape(permute(base_FixC,[1 3 2 4]),K*V,T*V);
Q= reshape(permute(base_Q,[1 3 2 4]),K*V,T*V);
%% Variables
purch= model.decision(M);
m= model.decision(M, T);
r= model.decision(M, T);
Iemp= model.decision(N, T);
v=  model.decision(K*V, V*T);
u=  model.decision(K*V, V*T);
x=  model.decision(K*V, V*T, 'B');
%
taw= model.decision(T);
obj3= model.decision(T);
obj4= model.decision(T);
obj5= model.decision(T);
obj7= model.decision(T);
%
obj1= model.decision(T-1);
obj2= model.decision(T-1);
obj6= model.decision(T-1);
%% Model Z
p= model.ambiguity;
p.suppset( A_D*delta <= b_D , A_P*pi <= b_P);
model.with(p);
vv=0;
uu=0;
for t=1:T
    r_sum=0;    
    %% Affine
    if t>=2
        inside1=[];
        for i=1:2*N
           inside1=[inside1, (i-1)*T+(1:t-1)]; 
        end
        r(:,t).affadapt(zeta(inside1));
        m(:,t).affadapt(zeta(inside1));
    end
    inside2=[];
    for i=1:N
        inside2=[inside2, (i-1)*T+(1:t)];
    end
    Iemp(:,t).affadapt(pi(inside2));
    v(:,(t-1)*V+1:(t-1)*V+V).affadapt(delta(inside2));
    u(:,(t-1)*V+1:(t-1)*V+V).affadapt(pi(inside2));
    %% End of Affine
    f=0;
    e=0;
    xx=0;
    for k=1:K
        f=f+ sum(v((k-1)*V+1:(k-1)*V+V, (t-1)*V+M+1:(t-1)*V+V),1)' - sum(v((k-1)*V+M+1:(k-1)*V+V, (t-1)*V+1:(t-1)*V+V),2) ; % related to Cons 19
        e=e+ sum(u((k-1)*V+M+1:(k-1)*V+V, (t-1)*V+1:(t-1)*V+V),2) - sum(u((k-1)*V+1:(k-1)*V+V, (t-1)*V+M+1:(t-1)*V+V),1)' ; % related to Cons 20
        %
        vv=vv +sum( v( (k-1)*V+1:(k-1)*V+M , (t-1)*V+M+1:(t-1)*V+V ) ,2);
        uu=uu +sum( u( (k-1)*V+M+1:(k-1)*V+V , (t-1)*V+1:(t-1)*V+M ) ,1);
        xx=xx +sum( x( (k-1)*V+M+1:(k-1)*V+V , (t-1)*V+1:(t-1)*V+V ), 2);        
        for i=M+1:V
           model.append( x((k-1)*V+i, (t-1)*V+i) == 0 );    % Routing from a retailer to itself is not allowed
        end    
        model.append( sum(sum(x((k-1)*V+1:(k-1)*V+M, (t-1)*V+M+1:(t-1)*V+V), 1), 2) <= 1 );        %Cons 16
        model.append( sum(sum(x( (k-1)*V+1:(k-1)*V+M, (t-1)*V+1:(t-1)*V+M ), 1), 2) == 0);       % Shipment among Breweries is not allowed
        model.append( sum(x((k-1)*V+1:(k-1)*V+V, (t-1)*V+1:(t-1)*V+V), 2) - sum(x((k-1)*V+1:(k-1)*V+V, (t-1)*V+1:(t-1)*V+V), 1)' == 0 );  %Cons 18
    end  
    if t==1        
        model.append( Iemp0(M+1:V,1) - e + pi(t:T:T*(N-1)+t) <= Iemp(1:N,t) );  %cons 9        
        model.append( r(:,1) <= alpha*Iemp0(1:M) );  %Cons 2
        model.append( m(:,t) <= purch );  %Cons 4
        Iemp_M= Iemp0(1:M,1) + uu' - ((1/alpha)*r(:,1));   % related to Cons 8
    else
        for s=1:t            
            r_sum=r_sum+ r(:,s);
        end
        model.append( Iemp(1:N,t-1) - e + pi(t:T:T*(N-1)+t) <= Iemp(1:N,t) );  %cons 9
        model.append( r(:,t) <= alpha*Iemp_prev_M );  %Cons 11
        model.append( m(:,t) <= purch-sum(m(:,1:t-1),2) );  %Cons 13
        Iemp_M= Iemp0(1:M,1) + uu' - ((1/alpha)*r_sum);   % related to Cons 8
    end
    model.append( f >= delta(t:T:T*(N-1)+t) );  %cons 10
    model.append( xx <= 1 );   %Cons 17
    model.append( Iemp_M <= Le(1:M) );  %Cons 15
    model.append( Iemp(1:N,t) <= Le(M+1:V) );  %Cons 15
    %
    Iful= Iful0 + sum(m(:,1:t),2) + sum(r(:,1:t),2) - vv;   %related to cons 7
    model.append( Iful <= Lf );  %Cons 14    
    model.append( Hf'*Iful <= taw(t) ); %Obj func
    %
    model.append( (He(1:M)'*Iemp_M) <= obj5(t) ); %Obj func
    model.append( (He(M+1:V)'*Iemp(1:N,t)) <= obj7(t) ); %Obj func    
    if t ~= T
        model.append( sum((Cins-((1-alpha)*Cdis))*Iemp_M) <= obj2(t) ); %Obj func
    end    
    model.append( m(:,t)+r(:,t) <= FilCap );  %Cons 3,12        
    model.append( Iful >= 0 ); %Cons 22
    model.append( f >= 0 ); %Cons 24
    model.append( e >= 0 ); %Cons 24
    model.append( vv >= 0 );
    model.append( uu >= 0 );
    model.append( xx >= 0 );
    model.append( Iemp_M >= 0 );    
    Iemp_prev_M= Iemp_M;        
end
model.append( v + u <= Q .* x );    %Cons 17
model.append( v >= 0 ); %Cons 25
model.append( u >= 0 ); %Cons 25
model.append( Iemp >= 0 ); %Cons 23
model.append( m >= 0 ); %Cons 22
model.append( r >= 0 ); %Cons 22
model.append( purch >= 0 ); %Cons 5
for t=2:T
    model.append( sum(Chyg.*r(:,t)) <= obj1(t-1) ); %Obj func
    model.append( (He(1:M)'*(purch-(sum(m(:,1:t),2)))) <= obj6(t-1) ); %Obj func
end
for t=1:T
    model.append( sum(sum(C(:,(t-1)*V+1:(t-1)*V+V).*(v(:,(t-1)*V+1:(t-1)*V+V)+u(:,(t-1)*V+1:(t-1)*V+V)))) <= obj3(t) ); %Obj func
    model.append( sum(sum(FixC(:,(t-1)*V+1:(t-1)*V+V).*x(:,(t-1)*V+1:(t-1)*V+V))) <= obj4(t) );      %Obj func
end
model.min( sum(Cpur.*purch)+ sum(Chyg.*r(:,1))+ sum((Cins-((1-alpha)*Cdis))*Iemp0(1:M))+ (He(1:M)'*(purch-m(:,1)))...
    + sum(obj1)...
    + sum(obj2)...
    + sum(obj3)...
    + sum(obj4)...
    + sum(obj5)...
    + sum(obj7)...
    + sum(obj6)...
    + sum(taw) )
%
model.solve
deter_time= model.Solution.time;
Obj= model.get;
solving_func_time= deter_time;
sol_status= model.Solution.status;

val.purch= purch.get;
val.m= m.get;
val.r= r.get;
% val.Iful= Iful;
% val.Iemp= Iemp.get;
% val.f= f;
% val.e= e;
% val.v= v.get;
% val.u= u.get;
val.x= x.get;
% val.delta= delta;
% val.pi= pi;
end