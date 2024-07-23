clc;
clear all;
X=linspace(0,1,200);
lambda=calc_lambda(10);
lambda=removeDuplicates(lambda);
lambda
N1=1;
N2=3;
theta1=zeros(size(X,2),1);
theta2=zeros(size(X,2),1);
theta3=zeros(size(X,2),1);
theta4=zeros(size(X,2),1);
theta5=zeros(size(X,2),1);
theta6=zeros(size(X,2),1);
for k=1:size(X,2)
    sum=0;
    for j=1:N1
        sum=sum+calc_theta(lambda(1,j),X(k),0.05);
    end
    theta1(k,1)=sum;
end
for k=1:size(X,2)
    sum=0;
    for j=1:N2
        sum=sum+calc_theta(lambda(1,j),X(k),0.05);
    end
    theta2(k,1)=sum;
end
for k=1:size(X,2)
    sum=0;
    for j=1:N1
        sum=sum+calc_theta(lambda(1,j),X(k),0.1);
    end
    theta3(k,1)=sum;
end
for k=1:size(X,2)
    sum=0;
    for j=1:N2
        sum=sum+calc_theta(lambda(1,j),X(k),0.1);
    end
    theta4(k,1)=sum;
end
for k=1:size(X,2)
    sum=0;
    for j=1:N1
        sum=sum+calc_theta(lambda(1,j),X(k),1);
    end
    theta5(k,1)=sum;
end
for k=1:size(X,2)
    sum=0;
    for j=1:N2
        sum=sum+calc_theta(lambda(1,j),X(k),1);
    end
    theta6(k,1)=sum;
end
tau=linspace(0,2,200);
X0=0;
Bi1=10;
Bi2=1;
theta7=zeros(size(tau,2),1);
theta8=zeros(size(tau,2),1);
for k=1:size(tau,2)
    sum=0;
    for j=1:N1
        sum=sum+calc_theta(lambda(1,j),X0,tau(1,k));
    end
    theta7(k,1)=sum;
end
for k=1:size(tau,2)
    sum=0;
    for j=1:N2
        sum=sum+calc_theta(lambda(1,j),X0,tau(1,k));
    end
    theta8(k,1)=sum;
end
theta9=zeros(size(tau,2),1);
theta10=zeros(size(tau,2),1);
lambda=calc_lambda(1);
lambda=removeDuplicates(lambda);
for k=1:size(tau,2)
    sum=0;
    for j=1:N1
        sum=sum+calc_theta(lambda(1,j),X0,tau(1,k));
    end
    theta9(k,1)=sum;
end
for k=1:size(tau,2)
    sum=0;
    for j=1:N2
        sum=sum+calc_theta(lambda(1,j),X0,tau(1,k));
    end
    theta10(k,1)=sum;
end
lambda=calc_lambda(0.1);
lambda=removeDuplicates(lambda);
theta_lumped=zeros(size(tau,2),1);
theta11=zeros(size(tau,2),1);
lambda
for k=1:size(tau,2)
    sum=0;
    for j=1:N2
        sum=sum+calc_theta(lambda(1,j),X0,tau(1,k));
    end
    theta11(k,1)=sum;
end
for i=1:size(tau,2)
theta_lumped(i,1)=calc_lumped(0.1,tau(1,i));
end
figure(1);
plot(X,theta1);
hold on; % Keep the current plot while adding new data
plot(X,theta2);
legend('N1=1','N2=3');
xlabel('X');
ylabel('Theta');
title('Plot of Theta for N1=1 and N2=3 at tau=0.05');
figure(2)
plot(X,theta3);
hold on;
plot(X,theta4);
legend('N1=1','N2=3'); % Add legend
xlabel('X');
ylabel('Theta');
title('Plot of Theta for tau=0.1');
figure(3)
plot(X,theta5);
hold on;
plot(X,theta6);
legend('N1=1','N2=3'); % Add legend
xlabel('X');
ylabel('Theta');
title('Plot of Theta for tau=1');
figure(4)
plot(tau,theta7);
hold on;
plot(tau,theta8);
xlabel('tau');
ylabel('Theta');
title('Plot of Theta vs tau for Bi=10');
legend('N=1','N=3');
figure(5)
plot(tau,theta9);
hold on;
plot(tau,theta10);
xlabel('tau');
ylabel('Theta');
title('Plot of Theta vs tau for Bi=1');
legend('N=1','N=3');
figure(6)
plot(tau,theta_lumped);
hold on;
plot(tau,theta11);
xlabel('tau');
ylabel('theta');
title('Plot of theta vs tau');
legend('lumped','non-lumped');

% Define the spatial mesh
x = linspace(0, 1, 100); % dimensionless spatial coordinate

% Define the time vector
t = linspace(0, 1, 100); % dimensionless time
N3=5;
theta12=zeros(size(x,2),1);
theta13=zeros(size(x,2),1);
theta14=zeros(size(x,2),1);
lambda=calc_lambda(10);
lambda=removeDuplicates(lambda);
for k=1:size(x,2)
    sum=0;
    for j=1:N3
        sum=sum+calc_theta(lambda(1,j),X(k),0.001);
    end
    theta12(k,1)=sum;
end
for k=1:size(x,2)
    sum=0;
    for j=1:N3
        sum=sum+calc_theta(lambda(1,j),X(k),0.01);
    end
    theta13(k,1)=sum;
end
for k=1:size(x,2)
    sum=0;
    for j=1:N3
        sum=sum+calc_theta(lambda(1,j),X(k),0.1);
    end
    theta14(k,1)=sum;
end
figure(7);
plot(x,theta12);
hold on;
plot(x,theta13);
hold on;
plot(x,theta14);
xlabel('X')
ylabel('theta')
legend('tau=0.001','tau=0.01','tau=0.1');
%% solved pdepe
clc;
clear all;
m = 0;
x = linspace(0, 1, 100);
t = linspace(0, 1, 1000);
sol = pdepe(m, @pdepefun, @initial_cond, @boundary_cond, x, t);
% Plot the solution
figure(1)
surf(x, t, sol(:,:,1));
xlabel('Distance (X)');
ylabel('Time (tau)');
zlabel('Solution (theta)');
title('Solution of the PDE');
theta1=sol(2,:,1);
theta2=sol(11,:,1);
theta3=sol(1000,:,1);
figure(2)
plot(x,theta1);
hold on;
plot(x,theta2);
hold on;
plot(x,theta3);
legend('tau=0.001','tau=0.01','tau=0.1');
xlabel('x');
ylabel('theta');
