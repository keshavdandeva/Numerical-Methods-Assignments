clear
close all
clc

A=[0 1 -1; 0 1 -2; 1200 -282 -62];
b=[0 ;0 ;-1];
c=[820 ;-296.5 ;-8];
ctrans= transpose(c);
x=@(t)1*(t>0); % case 1
x2=@(t)exp(-t)*(t>0); %case 2

v0=[0; 0; 0];

counter = 0;

for h=0.001:0.01:0.25   % For checking various values of integration step
    
    counter= counter +1;
    stepsize(counter)=h;
    
    t=0:h:5;
    %% ode45: Case 1
    dvdt=@(t,v) A*v+b*x(t);
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [to,vo] = ode45(dvdt,t,v0,options);
    
    votrans = transpose(vo);
    y = ctrans*votrans;
    
    
    %% ode45: Case 2
    
    dvdt2=@(t,v) A*v+b*x2(t);
    [to2,vo2] = ode45(dvdt2,t,v0,options);
    
    vo2trans = transpose(vo2);
    y2 = ctrans*vo2trans;
    
%         figure(1)       % This is the graph of ODE Solution
%         plot(t,y,'r-');
%         grid on
%         hold on
%         plot(t,y2,'g-');
%         title('ODE45 Method');
%         xlabel('t');
%         ylabel('y(t)');
%         legend('Case:1','Case:2');
%     
    %% Explicit Adams-Bashforth method: Case 1
    
    vabe(:,1)= v0; % K = 0
    vabe(:,2)= vabe(:,1)+h*(A*vabe(:,1)+b*x(t(1))); % K = 1
    
    for i=3:length(t)
        vabe(:,i)=vabe(:,i-1)+h*((3/2)*(A*vabe(:,i-1)+b*x(t(i-1)))+(-1/2)*(A*vabe(:,i-2)+b*x(t(i-2))));
    end
    y3 = ctrans*vabe;
    
    clear vabe
    
    abe1_Delta2(counter) = norm(y3-y,2);
    abe1_DeltaInf(counter)= norm(y3-y,Inf);
    
    %Task-3 Graphs
%     figure(4)                 
%     semilogy(stepsize,abe1_Delta2);
%     grid on
%     title('Adam Bashforth Method: Case 1');
%     xlabel('Step-size (h)');
%     ylabel('\delta_2(h)');
% 
%     figure(5)
%     semilogy(stepsize,abe1_DeltaInf);
%     grid on
%     title('Adam Bashforth Method: Case 1');
%     xlabel('Step-size (h)');
%     ylabel('\delta_{\infty}(h)');
    

    %% Explicit Adams-Bashforth method: Case 2
    
    vabe2(:,1)= v0; % K = 0
    vabe2(:,2)= vabe2(:,1)+h*(A*vabe2(:,1)+b*x2(t(1))); % K = 1
    
    for i=3:length(t)
        vabe2(:,i)=vabe2(:,i-1)+h*((3/2)*(A*vabe2(:,i-1)+b*x2(t(i-1)))+(-1/2)*(A*vabe2(:,i-2)+b*x2(t(i-2))));
    end
    
    y4 = ctrans*vabe2;
    clear vabe2
    
    abe2_Delta2(counter) = norm(y4-y2,2);
    abe2_DeltaInf(counter) = norm(y4-y2,Inf);
    
    %     figure(2)                 % This is the graph of Adam Bashforth Solution
    %     plot(t,y3,'-r');
    %     hold on
    %     plot(t,y4,'-g');
    %     grid on
    %     title('Explicit Adams-Bashforth Method');
    %     xlabel('t');
    %     ylabel('y(t)');
    %     legend('Case:1','Case:2');
    
    %Task 3 Graph
%     figure(6)
%     semilogy(stepsize,abe2_Delta2);
%     grid on
%     title('Adam Bashforth Method: Case 2');
%     xlabel('Step-size (h)');
%     ylabel('\delta_2(h)');
% 
%     figure(7)
%     semilogy(stepsize,abe2_DeltaInf);
%     grid on
%     title('Adam Bashforth Method: Case 2');
%     xlabel('Step-size (h)');
%     ylabel('\delta_{\infty}(h)');
    %% Implicit Gear Method: Case 1
    
    vig(:,1)= v0; %K=0
    vig(:,2)= (eye(3)-h*1*A)\(h*1*b*x(t(2)) + (1*vig(:,1))); %K=1
    vig(:,3)= (eye(3)-h*(2/3)*A)\(h*(2/3)*b*x(t(3)) + ((4/3)*vig(:,2)) + ((-1/3)*vig(:,1))); %K=2
    vig(:,4)= (eye(3)-h*(6/11)*A)\(h*(6/11)*b*x(t(4)) + ((18/11)*vig(:,3)) + ((-9/11)*vig(:,2)) + ((2/11)*vig(:,1))); %K=3
    vig(:,5)= (eye(3)-h*(12/25)*A)\(h*(12/25)*b*x(t(5)) + ((48/25)*vig(:,4)) + ((-36/25)*vig(:,3)) + ((16/25)*vig(:,2)) + ((-3/25)*vig(:,1))); %K=4
    
    beta0=60/137;
    alpha1=300/137;
    alpha2=-300/137;
    alpha3=200/137;
    alpha4=-75/137;
    alpha5=12/137;
    
    for i=6:length(t)
        vig(:,i)=(eye(3)-h*beta0*A)\(h*beta0*b*x(t(i)) + (alpha1*vig(:,i-1)) + (alpha2*vig(:,i-2)) + (alpha3*vig(:,i-3)) + (alpha4*vig(:,i-4)) + (alpha5*vig(:,i-5)));
    end
    
    y5 = ctrans*vig;
    clear vig
    ig1_Delta2(counter) = norm(y5-y,2);
    ig1_DeltaInf(counter) = norm(y5-y,Inf);
    
%     Task-3 Graphs
%     figure(8)
%     semilogy(stepsize,ig1_Delta2);
%     grid on
%     title('Gear Method: Case 1');
%     xlabel('Step-size (h)');
%     ylabel('\delta_2(h)');
%     
%     figure(9)
%     semilogy(stepsize,ig1_DeltaInf);
%     grid on
%     title('Gear Method: Case 1');
%     xlabel('Step-size (h)');
%     ylabel('\delta_{\infty}(h)');
%   
    
    %% Implicit Gear Method: Case 2
    
    vig2(:,1)= v0; %K=0
    vig2(:,2)= (eye(3)-h*1*A)\(h*1*b*x2(t(2)) + (1*vig2(:,1))); %K=1
    vig2(:,3)= (eye(3)-h*(2/3)*A)\(h*(2/3)*b*x2(t(3)) + ((4/3)*vig2(:,2)) + ((-1/3)*vig2(:,1))); %K=2
    vig2(:,4)= (eye(3)-h*(6/11)*A)\(h*(6/11)*b*x2(t(4)) + ((18/11)*vig2(:,3)) + ((-9/11)*vig2(:,2)) + ((2/11)*vig2(:,1))); %K=3
    vig2(:,5)= (eye(3)-h*(12/25)*A)\(h*(12/25)*b*x2(t(5)) + ((48/25)*vig2(:,4)) + ((-36/25)*vig2(:,3)) + ((16/25)*vig2(:,2)) + ((-3/25)*vig2(:,1))); %K=4
    
    for i=6:length(t)
        vig2(:,i)=(eye(3)-h*beta0*A)\(h*beta0*b*x2(t(i)) + (alpha1*vig2(:,i-1)) + (alpha2*vig2(:,i-2)) + (alpha3*vig2(:,i-3)) + (alpha4*vig2(:,i-4)) + (alpha5*vig2(:,i-5)));
    end
    
    y6 = ctrans*vig2;
    clear vig2
    
    ig2_Delta2(counter) = norm(y6-y2,2);
    ig2_DeltaInf(counter) = norm(y6-y2,Inf);
    
    %     figure(3)                  % This is the graph of Implicit Gear Method Solution
    %     plot(t,y5,'r-');
    %     hold on
    %     plot(t,y6,'g-');
    %     grid on
    %     title('Implicit Gear Method');
    %     xlabel('t');
    %     ylabel('y(t)');
    %     legend('Case:1','Case:2');
    
    %     Task-3 Graphs
%     figure(10)
%     semilogy(stepsize,ig2_Delta2);
%     grid on
%     title('Gear Method: Case 2');
%     xlabel('Step-size (h)');
%     ylabel('\delta_2(h)');
%     
%     figure(11)
%     semilogy(stepsize,ig2_DeltaInf);
%     grid on
%     title('Gear Method: Case 2');
%     xlabel('Step-size (h)');
%     ylabel('\delta_{\infty}(h)');
    
    
end


