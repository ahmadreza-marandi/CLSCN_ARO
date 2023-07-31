function [newA, newb, solving_func_time] = Alg7_func (A_ini ,b_ini, D, incornerL, incornerU, prelambda, lambda)
%Partitioning
%% Parameters
b_ini= round(b_ini, 10);
A_row_num=size(A_ini,1);
left=1;
right=2;
newA=double.empty(1,0);
newb=double.empty(1,0);
teek0=tic;
%% Main
xc=(incornerU + incornerL)/2;
difc=(incornerU - incornerL)/2;
for k=1:D
    for i=1:A_row_num
        if A_ini(i,k) < 0
            prelambda(2*k-1)= (b_ini(i)-(A_ini(i,:)*xc))/A_ini(i,k);
            if prelambda(2*k-1) > lambda(2*k-1)
                lambda(2*k-1)= prelambda(2*k-1);
            end
        elseif A_ini(i,k) > 0
            prelambda(2*k)= (b_ini(i)-(A_ini(i,:)*xc))/A_ini(i,k);
            if prelambda(2*k) < lambda(2*k)
                lambda(2*k)= prelambda(2*k);
            end
        end
    end
    lambda(2*k-1)= lambda(2*k-1) + difc(k);       %% 21
    lambda(2*k)= lambda(2*k) - difc(k);       %% 21        
    if lambda(2*k-1) < -0.1   % -100: we don't want this part because of first polytope shape
        newA(left).part= [A_ini; zeros((k-1)*2+1,D) ];
        newb(left).part= [b_ini; zeros((k-1)*2+1,1) ];        
        for j =1:k-1  % for k >1      
            pikRow= A_row_num+((j-1)*2+1);
            newA(left).part(pikRow, j)= -1;
            newA(left).part(pikRow+1, j)= 1;
            newb(left).part(pikRow:pikRow+1, 1)= [ -incornerL(j); incornerU(j) ];
        end
        newA(left).part(end, k)= 1;
        newb(left).part(end, 1)= incornerL(k);
        left=right+1;
    else
        right=right-1;
        left=right+1;
    end        
    if lambda(2*k) > 0.01        
        newA(right).part= [A_ini; zeros((k-1)*2+1,D) ];
        newb(right).part= [b_ini; zeros((k-1)*2+1,1) ];            
        for j =1:k-1  % for k >1      
            pikRow= A_row_num+((j-1)*2+1);
            newA(right).part(pikRow, j)= -1;
            newA(right).part(pikRow+1, j)= 1;
            newb(right).part(pikRow:pikRow+1, 1)= [ -incornerL(j); incornerU(j) ];
        end
        newA(right).part(end, k)= -1;
        newb(right).part(end, 1)= -incornerU(k);
        right=left+1;
    else
        left=left-1;
        right=left+1;
    end
end
solving_func_time=toc(teek0);
end