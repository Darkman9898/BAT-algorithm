%Author: [Thanh-Tung Nguyen]
%Email: tungnguyen98ac@gmail.com
%Workplace: Ho Chi Minh City Bach Nghe College, Ho Chi Minh City (HCMC), Vietnam
%Date: June 08, 2025
%Description: MATLAB code for optimizing LQR controller parameters using the Bat algorithm.

clc;
clear all;
close all;

%A=6x6
%B=6x1

para = [20 100 0.9 0.7]; 
n = para(1);              
N_gen = para(2);         
AA = para(3);             
r = para(4);             


Qmin = 0;                 
Qmax = 2;                 


d = 7;                   

Lb = [0 0 0 0 0 0 1];     
Ub = [1000000, 100000000, 100000000, 100000000, 100000000, 100000000, 10];


Q = zeros(n, 1);         
v = zeros(n, d);          
Sol = zeros(n, d);        
Fitness = zeros(n, 1);    


for i = 1:n
    Sol(i, :) = Lb + (Ub - Lb) .* rand(1, d); 
end


Tmax = 50;               
tck = 0.005;             


best = Sol(1, :);        
fmin = Inf;               

for i = 1:n
    
    if any(Sol(i, 1:7) <= 0)  
        Fitness(i) = Inf;      
        continue;
    end

    
    q = diag([Sol(i, 1), Sol(i, 2), Sol(i, 3), Sol(i, 4), Sol(i, 5), Sol(i, 6)]);
    r = Sol(i, 7);
    assignin('base', 'r', r);
    assignin('base', 'q', q);
    [K, ~, ~] = lqr(A, B, q, r);
    assignin('base', 'K', K);


    model_name = 'RDPIP_LQR_Simulink';
    try
        simOut = sim(model_name, 'StopTime', '50', ...
                     'SaveOutput', 'on', 'SrcWorkspace', 'base');
        
        
        e1 = simOut.get('e1');
        if length(e1) > Tmax / tck
            e2 = simOut.get('e2');
            e3 = simOut.get('e3');
            Fitness(i) = sum(e1) + 1 * sum(e2) + 1 * sum(e3) + 1;
        else
            Fitness(i) = Inf;
        end
    catch
        Fitness(i) = Inf;
    end

    if Fitness(i) < fmin
        best = Sol(i, :);
        fmin = Fitness(i);
    end
end

for t = 1:N_gen
    for i = 1:n
        Q(i) = Qmin + (Qmax - Qmin) * rand;
        v(i, :) = v(i, :) + (Sol(i, :) - best) * Q(i);
        S(i, :) = Sol(i, :) + v(i, :);

        S(i, :) = max(S(i, :), Lb);
        S(i, :) = min(S(i, :), Ub);

        if rand > r
            S(i, :) = best + 0.001 * randn(1, d);
        end

        if any(S(i, 1:7) <= 0)
            Fnew = Inf;
        else
            q = diag([S(i, 1), S(i, 2), S(i, 3), S(i, 4), S(i, 5), S(i, 6)]);
            r = S(i, 7);
            assignin('base', 'r', r);
            assignin('base', 'q', q);
            [K, ~, ~] = lqr(A, B, q, r);
            assignin('base', 'K', K);

            try
                simOut = sim(model_name, 'StopTime', '50', ...
                             'SaveOutput', 'on', 'SrcWorkspace', 'base');
                
                e1 = simOut.get('e1');
                if length(e1) > Tmax / tck
                    e2 = simOut.get('e2');
                    e3 = simOut.get('e3');
                    Fnew = sum(e1) + 1 * sum(e2) + 1 * sum(e3) + 1;
                else
                    Fnew = Inf;
                end
            catch
                Fnew = Inf;
            end
        end

        if (Fnew <= Fitness(i)) && (rand < AA)
            Sol(i, :) = S(i, :);
            Fitness(i) = Fnew;
            AA = 0.9 * AA;
            r = r * (1 - exp(-0.1 * t));
        end

        if Fnew <= fmin
            best = S(i, :);
            fmin = Fnew;
        end
    end
end

disp(['Nghiệm tốt nhất: ', num2str(best)]);
disp(['Giá trị hàm mục tiêu tốt nhất: ', num2str(fmin)]);

