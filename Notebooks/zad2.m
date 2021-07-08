close all
clear all
clc

%% Parametri genetskog algoritma.

N = 1500; L = 3*24; pc = 0.9;
pm = 0.001; G = 0.8; C = 151;

%% Inicijalizacija populacije

gen_num = 0;

for i=1:N
    gen(i, 1:L) = dec2bin(round(rand*(2^L-1)), L);
end

[fitness_init, f1, f2, f3] = evaluiraj_populaciju_h(gen);
[max_fitness, max_ind] = max(fitness_init);
display(['Prosecna vrednost funkcije greske: ', num2str(mean(fitness_init))]);
display(['Maksimalna vrednost funkcije greske:', num2str(max_fitness)]);

[max_fitness, max_ind] = max(f1);
display(['Prosecna vrednost funkcije greske za prvu jednacinu: ', num2str(mean(f1))]);
display(['Maksimalna vrednost funkcije greske za prvu jednacinu:', num2str(max_fitness)]);

[max_fitness, max_ind] = max(f2);
display(['Prosecna vrednost funkcije greske za drugu jednacinu: ', num2str(mean(f2))]);
display(['Maksimalna vrednost funkcije greske za drugu jednacinu:', num2str(max_fitness)]);

[max_fitness, max_ind] = max(f3);
display(['Prosecna vrednost funkcije greske za trecu jednacinu: ', num2str(mean(f3))]);
display(['Maksimalna vrednost funkcije greske za trecu jedncinu:', num2str(max_fitness)]);

% x1 = bin2dec(gen(max_ind, 1:L/3)) * 2*pi / (2^(L/3)-1);
% y1 = bin2dec(gen(max_ind, L/3+1:2*L/3)) * 2*pi / (2^(L/3)-1);
% z1 = bin2dec(gen(max_ind, 2*L/3:end)) * 2*pi / (2^(L/3)-1);
% 
% display(['Najbolja resenja X = ', x1, 'Y = ', y1, 'Z = ', z1]);

%% Genetski algortam.
uslov = 0;

while uslov == 0
    
    fitness = evaluiraj_populaciju(gen);

    [fitness_max, ind_max] = max(fitness);

    display(['Najmanja greska za generaciju ', num2str(gen_num), ' je ', num2str(151 - fitness_max)]);
    
    fitness_sortirano = sort(fitness);
    if (gen_num > 350) || (abs(fitness_sortirano(N) - fitness_sortirano(N-10)) < 0.0001)
        
        [fitness_max, ind_max] = max(fitness);
        x1 = bin2dec(gen(max_ind, 1:L/3)) * 4 / (2^(L/3)-1);
        y1 = bin2dec(gen(max_ind, L/3+1:2*L/3)) * 4 / (2^(L/3)-1);
        z1 = bin2dec(gen(max_ind, 2*L/3:end)) * 4 / (2^(L/3)-1);
        
        [f1, f2, f3] = resi_jednacine(x1, y1, z1);
        eval = evaluiraj_populaciju(gen(max_ind, 1:L));
        display(num2str(C - eval));
        display(['Konacne vrednosti: x = ', num2str(x1), ' y = ', num2str(y1), ' z = ', num2str(z1)]);
        display(['Greske za jednacine : f1 = ', num2str(f1), ' f2 = ', num2str(f2), ' f3 = ', num2str(f3)]);
        break;
    end

    %% Formiranje nove generacije.

    N_ukrstanje = round(G*N);
    N_ukrstanje = N_ukrstanje + mod(N_ukrstanje, 2);

    N_reprodukcija = N - N_ukrstanje;

    % Reprodukcija.
    [s, ind_s] = sort(fitness, 'ascend');
    nova_gen(1:3, 1:L) = gen(ind_s(1:3), 1:L);
    nova_gen(4:N_reprodukcija, 1:L) = selektuj_tocak_ruleta(gen, fitness, N_reprodukcija-3);

    % Ukrstanje.
    cur_ind = N_reprodukcija + 1;
    for i=1:N_ukrstanje/2
        roditelji = selektuj_tocak_ruleta(gen, fitness, 2);
        nova_gen(cur_ind:cur_ind+1, 1:L) = ukrstanje(roditelji, pc);
        cur_ind = cur_ind + 2;
    end

    % Mutacija.
    N_mutiranih = min(1, round(N*L*pm));
    for i=1:N_mutiranih
        Nm = ceil(rand*N);
        Lm = ceil(rand*L);
        nova_gen(Nm, Lm) = num2str(1 - str2num(nova_gen(Nm, Lm)));
    end
    
    gen = nova_gen;
    gen_num = gen_num + 1;
end