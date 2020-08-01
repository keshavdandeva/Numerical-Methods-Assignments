clear
close all
clc
format long

% TASK:1

y =@(x) (6/5).^x + x + sin(x.^(1/2)) - 33/4 ;
x0 = [0 10];
Exact_Solution = fzero(y,x0);

figure(1)
x = linspace(0,10,100);
y = (6/5).^x + x + sin(x.^(1/2)) - 33/4 ;
plot(x,y);
grid on;
hold on;
plot(Exact_Solution,0,'-o');
text(Exact_Solution,0,'$\leftarrow \dot{x}$', 'Interpreter','latex');
xlabel('x');
ylabel('f(x)');
title('Plot of the Function');

% Task:2


for i = (1:14)
    delta = 10.^-(i+2);
    Bisection_epsilon(i) = delta;

    f =@(x) (6/5).^x + x + sin(x.^(1/2)) - 33/4 ;
    a = 10^(-6);
    b = 10;
    xi = (a + b)./2;
    counter = 0;

    if f(a)*f(b)>0
        disp('Error! Values provided does not satisfy Bisection method')
    else
        while abs((b-a)./2) > delta
            if f(a)*f(xi)<0
                b = xi;
            else
                a = xi;
            end
            xi = (a + b)/2;
            counter = counter + 1;

        end
        Bisection_StoppingValue(i) = abs((b-a)./2);
        Bisection_Roots(i) = xi;
        Bisection_AbsoluteError(i) = abs(Exact_Solution - xi);
        Bisection_Iterations(i) = counter;
    end
end

figure(2);
semilogx(Bisection_epsilon,Bisection_Iterations,'--x');
grid on
legend('Bisection')
xlabel('Epsilon');
ylabel('Number of Iterations');
title('I vs Epsilon');

figure(3);
loglog(Bisection_epsilon,Bisection_AbsoluteError);
hold on
loglog(Bisection_epsilon,Bisection_StoppingValue);
hold off
legend('Absolute Error','Stopping Criteria')
xlabel('Epsilon');
title('Bisection Method');

% Task:3


for i=1:14
    delta2 = 10.^-(i+2);
    RFM_epsilon(i) = delta2;
    f =@(x) (6/5).^x + x + sin(x.^(1/2)) - 33/4 ;
    x0 = 6;
    f0 = f(x0);
    x1 = 0;
    f1 = f(x1);
    count = 0;
    while count < 100
        x2 = x1 -((x1-x0)./(f1-f0)).*f1;

        if abs(x2-x1) < delta2
            break;
        end

        x1 = x2;
        f1 = f(x1);
        count = count + 1;
    end
    RFM_StoppingValue(i) = abs(x2-x1);
    RFM_roots(i) = x2;
    RFM_AbsoluteError(i) = abs(Exact_Solution - x2);
    RFM_Iterations(i) = count;
end

figure(2);
semilogx(Bisection_epsilon,Bisection_Iterations,'--x');
grid on
hold on
semilogx(RFM_epsilon,RFM_Iterations,'--s');
hold off
legend('Bisection','RegulaFalsi')
xlabel('Epsilon');
ylabel('Number of Iterations');
title('I vs Epsilon');

figure(4);
loglog(RFM_epsilon,RFM_AbsoluteError);
hold on
loglog(RFM_epsilon,RFM_StoppingValue);
hold off
legend('Absolute Error','Stopping Criteria')
xlabel('Epsilon');
title('Regula Falsi Method');

% Task:4


for i = 1:14
    delta3 = 10^-(i+2);
    Secant_epsilon(i) = delta3;
    
    f =@(x) (6/5).^x + x + sin(x.^(1/2)) - 33/4 ;
    x0 = 0;
    x1 = 10;
    
    count = 0;
    while true
        x2 = x1 - (((x1-x0).*f(x1))./(f(x1)-f(x0)));
        
        if abs(x1-x0) < delta3
            break;
        end
        count = count + 1;
        x0 = x1;
        x1 = x2;
        
    end
    Secant_Iterations(i) = count;
    Secant_StoppingValue(i) = abs(x1-x0);
    Secant_roots(i) = x2;
    Secant_AbsoluteError(i) = abs(Exact_Solution - x2);
end


figure(2);
semilogx(Bisection_epsilon,Bisection_Iterations,'--x');
grid on
hold on
semilogx(RFM_epsilon,RFM_Iterations,'--s');
semilogx(Secant_epsilon,Secant_Iterations,'--o');
hold off
legend('Bisection','RegulaFalsi','Secant')
xlabel('Epsilon');
ylabel('Number of Iterations');
title('I vs Epsilon');

figure(5);
loglog(Secant_epsilon,Secant_AbsoluteError);
hold on
loglog(Secant_epsilon,Secant_StoppingValue);
hold off
legend('Absolute Error','Stopping Criteria')
xlabel('Epsilon');
title('Secant Method');

% Task:5

f =@(x) (6/5).^x + x + sin(x.^(1/2)) - 33/4 ;
df =@(x) ((log(6).*6.^x)-(log(5).*6.^x))./5.^x + (cos(x.^(1/2)))./(2.*x.^(1/2)) + 1;


for i = 1:14
    delta4 = 10^-(i+2);
    Newton_epsilon(i) = delta4;
    x0 = 6;
    
    count = 0;
    while count < 100
        x1 = x0 - (f(x0)./df(x0));
        if abs(x1-x0) < delta4
            break;
        end
        temp_x0 = x0;
        x0 = x1;
        
        count = count + 1;
    end
    Newton_Iterations(i) = count;
    Newton_StoppingValue(i) = abs(x1-temp_x0);
    Newton_roots(i) = x1;
    Newton_AbsoluteError(i) = abs(Exact_Solution - x1);
end

figure(6);
loglog(Newton_epsilon,Newton_AbsoluteError);
hold on
loglog(Newton_epsilon,Newton_StoppingValue);
hold off
legend('Absolute Error','Stopping Criteria')
xlabel('Epsilon');
title('Newton Method');

figure(2);
semilogx(Bisection_epsilon,Bisection_Iterations,'--x');
grid on
hold on
semilogx(RFM_epsilon,RFM_Iterations,'--s');
semilogx(Secant_epsilon,Secant_Iterations,'--o');
semilogx(Newton_epsilon,Newton_Iterations,'--*');
hold off
legend('Bisection','RegulaFalsi','Secant','Newton')
xlabel('Epsilon');
ylabel('Number of Iterations');
title('I vs Epsilon');




