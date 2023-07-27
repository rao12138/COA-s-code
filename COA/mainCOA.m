%  Crayfish Optimization Algorithm(COA)
%
%  Source codes demo version 1.0                                                                      
%                                                                                                     
%  The 11th Gen Intel(R) Core(TM) i7-11700 processor with the primary frequency of 2.50GHz, 16GB memory, and the operating system of 64-bit windows 11 using matlab2021a.                                                                
%                                                                                                     
%  Author and programmer: Heming Jia,Honghua Rao,Changsheng Wen,Seyedali Mirjalili                                                                          
%         e-Mail: jiaheminglucky99@126.com;rao12138@163.com                                                                                                                                                                                                                                                


clear all 
close all
clc

N=100;  %Number of search agents
F_name='F3';     %Name of the test function
T=200;           %Maximum number of iterations

    
[lb,ub,dim,fobj]=Get_F(F_name); %Get details of the benchmark functions
[best_fun,best_position,cuve_f,global_Cov]=COA(N,T,lb,ub,dim,fobj); 


figure('Position',[454   445   694   297]);
subplot(1,3,1);
func_plot(F_name);     % Function plot
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([F_name,'( x_1 , x_2 )'])
subplot(1,3,2);       % Convergence plot
semilogy(cuve_f,'LineWidth',3)
subplot(1,3,3);       % Convergence plot
semilogy(global_Cov,'LineWidth',3)
xlabel('Iteration#');
ylabel('Best fitness so far');
legend('COA');



display(['The best-obtained solution by COA is : ', num2str(best_position)]);  
display(['The best optimal value of the objective funciton found by COA is : ', num2str(best_fun)]);  