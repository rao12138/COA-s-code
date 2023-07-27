%  Crayfish Optimization Algorithm(COA)
%
%  Source codes demo version 1.0                                                                      
%                                                                                                     
%  The 11th Gen Intel(R) Core(TM) i7-11700 processor with the primary frequency of 2.50GHz, 16GB memory, and the operating system of 64-bit windows 11 using matlab2021a.                                                                
%                                                                                                     
%  Author and programmer: Heming Jia,Raohonghua,Wen changsheng,Seyedali Mirjalili                                                                          
%         e-Mail: jiaheminglucky99@126.com;rao12138@163.com                                                           
%                                   
%                                                                                                                                                   
%_______________________________________________________________________________________________
% You can simply define your cost function in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% T = maximum number of iterations
% N = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single numbers
function [best_fun,best_position,cuve_f,global_Cov]  =COA(N,T,lb,ub,dim,fobj)
%% Define Parameters
cuve_f=zeros(1,T); 
X=initialization(N,dim,ub,lb); %Initialize population
global_Cov = zeros(1,T);
Best_fitness = inf;
best_position = zeros(1,dim);
fitness_f = zeros(1,N);
for i=1:N
   fitness_f(i) =  fobj(X(i,:)); %Calculate the fitness value of the function
   if fitness_f(i)<Best_fitness
       Best_fitness = fitness_f(i);
       best_position = X(i,:);
   end
end
global_position = best_position; 
global_fitness = Best_fitness;
cuve_f(1)=Best_fitness;
t=1; 
while(t<=T)
    C = 2-(t/T); %Eq.(7)
    temp = rand*15+20; %Eq.(3)
    xf = (best_position+global_position)/2; % Eq.(5)
    Xfood = best_position;
    for i = 1:N
        if temp>30
            %% summer resort stage
            if rand<0.5
                Xnew(i,:) = X(i,:)+C*rand(1,dim).*(xf-X(i,:)); %Eq.(6)
            else
            %% competition stage
                for j = 1:dim
                    z = round(rand*(N-1))+1;  %Eq.(9)
                    Xnew(i,j) = X(i,j)-X(z,j)+xf(j);  %Eq.(8)
                end
            end
        else
            %% foraging stage
            P = 3*rand*fitness_f(i)/fobj(Xfood); %Eq.(4)
            if P>2   % The food is too big
                 Xfood = exp(-1/P).*Xfood;   %Eq.(12)
                for j = 1:dim
                    Xnew(i,j) = X(i,j)+cos(2*pi*rand)*Xfood(j)*p_obj(temp)-sin(2*pi*rand)*Xfood(j)*p_obj(temp); %Eq.(13)
                end
            else
                Xnew(i,:) = (X(i,:)-Xfood)*p_obj(temp)+p_obj(temp).*rand(1,dim).*X(i,:); %Eq.(14)
            end
        end
    end
    %% boundary conditions
    for i=1:N
        for j =1:dim
            if length(ub)==1
                Xnew(i,j) = min(ub,Xnew(i,j));
                Xnew(i,j) = max(lb,Xnew(i,j));
            else
                Xnew(i,j) = min(ub(j),Xnew(i,j));
                Xnew(i,j) = max(lb(j),Xnew(i,j));
            end
        end
    end
   
    global_position = Xnew(1,:);
    global_fitness = fobj(global_position);
 
    for i =1:N
         %% Obtain the optimal solution for the updated population
        new_fitness = fobj(Xnew(i,:));
        if new_fitness<global_fitness
                 global_fitness = new_fitness;
                 global_position = Xnew(i,:);
        end
        %% Update the population to a new location
        if new_fitness<fitness_f(i)
             fitness_f(i) = new_fitness;
             X(i,:) = Xnew(i,:);
             if fitness_f(i)<Best_fitness
                 Best_fitness=fitness_f(i);
                 best_position = X(i,:);
             end
        end
    end
    global_Cov(t) = global_fitness;
    cuve_f(t) = Best_fitness;
    t=t+1;
    if mod(t,50)==0
      disp("COA"+"iter"+num2str(t)+": "+Best_fitness); 
   end
end
 best_fun = Best_fitness;
end
function y = p_obj(x)   %Eq.(4)
    y = 0.2*(1/(sqrt(2*pi)*3))*exp(-(x-25).^2/(2*3.^2));
end