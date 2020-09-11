# fea-1D-bar
Code written in MATLAB to solve FEA 1D Linear bar element problems 
Created as a Mini Project with fellow contributors: 

* Guru Prasaath P 
* David Smith Sundarsingh 
* Gopi Karthick R 

### Script file: 
Linear_bar_selfWeight_final2.m is the script file which has the finite element code 
### Code: 
```clc;
%% Inputs
g = 9.81;
E = input('Enter youngs modulus in N/mm^2 (constant)');
D = input('Enter density in kg/m^3 (constant)');
D = D * 10^-9;
A = input('Enter Area of bar elements in mm^2(constant) ');
 
%% length inputs
m = input('Number of elements ');
L = zeros(m,1);
for i=1:m
    fprintf('Enter length of element %d in mm ',i);
    L(i,1)= input('');
end
%% global stiffness matrix
KGlobal=zeros(m+1,m+1);
for i = 1:m
    KLocal = ((A*E)/L(i,1)).*[1 -1;-1 1];
    KGlobal(i,i) = KGlobal(i,i)+KLocal(1,1);
    KGlobal(i,i+1) = KGlobal(i,i+1)+KLocal(1,2);
    KGlobal(i+1,i) = KGlobal(i+1,i)+KLocal(2,1);
    KGlobal(i+1,i+1) = KGlobal(i+1,i+1)+KLocal(2,2);
end
disp('Global stifness matrix, G = ');
disp(KGlobal);
%% Force matrix
FGlobal = zeros(m+1,1);
f = zeros(m+1,1);
for i = 1:m+1
    fprintf('Enter Force at node (%d) in N ',i);
    f(i,1) = input('');
end
for i = 1:m
    FLocal = (D*g*0.5*A*L(i,1)).*[1;1];
    FGlobal(i,1)= FGlobal(i,1)+ FLocal(1,1);
    FGlobal(i+1,1) = FGlobal(i+1,1) + FLocal(2,1);
end
FGlobal = -1.*FGlobal + f;
disp('Global Force matrix, F = ');
disp(FGlobal);
 
%% Displacements matrix
count = 0;
fprintf('Enter boundary conditions for nodal displacements:');
u = ones(m+1,1);
fprintf('\nWhich nodes are fixed? input row array ');
x = input('')';
xL = length(x);
%validity check
prompt1 = ('Unexpected no. of boundary conditions... please rerun code section ');
prompt2 = ('Unexpected boundary condition.. resuming anyway...');
if xL>=m+1
    disp(prompt1);
    return;
end
for j=1:xL
    if x(j)<=m+1
        u(x(j)) = 0;
    else
        if count==0;
            u(1)=0;
            disp(prompt2);
            break;;
        else
            disp(prompt2);
            continue;
        end
    end
    count = count+1;
end
 
for i=1:m+1
    if u(i)==0
        KGlobal(i,1:m+1) = 0;
        KGlobal(1:m+1,i) = 0;
        FGlobal(i,1) = 0;
    end
end
fprintf('\nAfter applying boundary conditions\n');
fprintf('K = \n');
disp(KGlobal);
fprintf('F = \n');
disp(FGlobal);
X = pinv(KGlobal)*FGlobal;
for i=1:m+1
    fprintf('U%d = %.3fmm\n',i,X(i));
end```




