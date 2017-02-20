clear
clc

%We want to find a model function m(X)=f+g^T*X+1/2*X^T*A*X
%Simply we want to find the i)scalar f, ii)vector g, and iii) Matrix A such
%that for q points y_i, m(y_i)=f(y_i); In this probelm y_i points are given
%and the function values f(y_i) at those points are given as well.

N=2;            %This is the dimension of the input array
q=(N+1)*(N+2)/2;%This is the number of constraints (equations)
Y=zeros(q,N);   %This is the matrix of the inputs, q points each N variables

Y(1,:)=[0 0];
Y(2,:)=[1 0];
Y(3,:)=[2 0];
Y(4,:)=[1 1];
Y(5,:)=[0 2];
Y(6,:)=[0 1];
f(1)=1;
f(2)=2.0084;
f(3)=7.0091;
f(4)=1.0168;
f(5)=-0.9909;
f(6)=-.9916;

% we want to build an equation Px=b where x is a vector of unknowns which:
% i) has the N*(N-1)/2 elements of entities of the symmteric matrix A
% ii)has the N unknowns values of g
% iii)has the scalar unknown value of f
% In total we have N(N-1)/2+N+1=(N+1)(N+2)/2=q unknonws.
% We are given q points y_1, y_2, ..., y_q and their functions values
% So have to solve a Linear System of Equation with q unknows P*x=b
% Need to find rows of matrix P such that P(i,:)'*x=b
% Here x is a vector arranged as
% x=[f,g_1,g_2,...,g_n,A_11,A12,...,A1N, A22,A23,...,A2N, A33,...,  ,ANN]
% For each rows of P,
% P(i,1)=1 --> f coefficient
% P(i,2:1+N)= [yi(1) yi(2), ..., yi(N)] --> g'*yi
% The following is based on X'*A*X
% P(i,2+N+1:2+N+1+N)=[yi(1)^2, 2*yi(1)*yi(2), ...., 2*yi(1)*yi(N)];
% P(i,4+2*N:4+2*N+(N-1))=[yi(2)^2, 2*yi(2)*yi(3), ...., 2*yi(2)*yi(N)];
% So on until the last one P(i,q)=yi(N)^2
% Will do it in a double for loop

P=zeros(q,q);
b=zeros(q,1);

for i=1:q %on coloumns
    P(i,1)=1;
    P(i,2:1+N)=Y(i,:);
    IDX_OFFSET=2+N;
    b(i,1)=f(i);
    % The double for loop to append A_ii coefficient to P(i,...)
    for j=1:N
        for k=j:N
            if(k~=j)
                value=Y(i,j)*Y(i,k);
            else
                value=1/2*Y(i,j)^2;
            end
            P(i,IDX_OFFSET)=value;
            IDX_OFFSET=IDX_OFFSET+1;
        end
    end
end

X=inv(P)*b



% Just to give you an idea of a simple 2D case above 
y1=[0 0]';
y2=[1 0]';
y3=[2 0]';
y4=[1 1]';
y5=[0 2]';
y6=[0 1]';
f1=1;
f2=2.0084;
f3=7.0091;
f4=1.0168;
f5=-0.9909;
f6=-.9916;

A=[1 y1(1) y1(2) 0.5*y1(1)^2  y1(1)*y1(2) 0.5*y1(2)^2;
    1 y2(1) y2(2) 0.5*y2(1)^2  y2(1)*y2(2) 0.5*y2(2)^2;
    1 y3(1) y3(2) 0.5*y3(1)^2  y3(1)*y3(2) 0.5*y3(2)^2;
    1 y4(1) y4(2) 0.5*y4(1)^2  y4(1)*y4(2) 0.5*y4(2)^2;
    1 y5(1) y5(2) 0.5*y5(1)^2  y5(1)*y5(2) 0.5*y5(2)^2;
    1 y6(1) y6(2) 0.5*y6(1)^2  y6(1)*y6(2) 0.5*y6(2)^2];

inv(A)*[f1 f2 f3 f4 f5 f6]'