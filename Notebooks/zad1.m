close all
clear all
clc

%% Generisanje i prikaz klasa.

N = 500;

X = zeros(3 * N, 2);

% Generisanje prve klase.
for i=1:N
    if rand < 0.5
        % Generisemo prvi kvadrat.
        X(i, 1) = rand;
        X(i, 2) = rand;
    else 
        % Generisemo drugi kvadrat.
        X(i, 1) = rand - 1;
        X(i, 2) = rand - 1;
    end
end

% Generisanje druge klase.
for i=N+1:2*N
    r = rand * 0.5;
    theta = rand * 2 * pi;
    X(i, 1) = r*cos(theta) - 0.5;
    X(i, 2) = r*sin(theta) + 0.5;
end

% Generisanje trece klase.
for i=2*N+1:3*N
    r = rand * 0.5;
    theta = rand * 2 * pi;
    X(i, 1) = r*cos(theta) + 0.5;
    X(i, 2) = r*sin(theta) - 0.5;
end

%% Prikaz generisanih podataka.

figure,
plot(X(1:N, 1), X(1:N, 2), 'r*');
hold on
plot(X(N+1:2*N, 1), X(N+1:2*N, 2), 'b*');
plot(X(2*N+1:3*N, 1), X(2*N+1:3*N, 2), 'g*');
hold off

%% Formiranje i obucavanje neuralne mreze.

X = X';

% One hot encoding - izlaz.
T = [ones(1, N), zeros(1, N), zeros(1, N);
     zeros(1, N), ones(1, N), zeros(1, N);
     zeros(1, N), zeros(1, N), ones(1, N)];
 
% [3 3] - podobucena
% [10 10] - mozda i preobucena

mreza = feedforwardnet([3 3]);
mreza.layers{1}.transferFcn = 'tansig';
mreza.layers{2}.transferFcn = 'tansig';
mreza.layers{3}.transferFcn = 'softmax';

mreza.trainParam.epochs = 1000;
mreza.trainParam.goal = 0.00001;

[mreza tr] = train(mreza, X, T);
err = min(tr.tperf)

%% Evaluacija i prikaz rezultata obucene mreze.

Y = sim(mreza, X);

k1 = 0; k2 = 0; k3 = 0;

MK = zeros(3, 3);

for i=1:3*N
    
    [max_y, ind] = max(Y(:, i));
    
    if ind == 1
        k1 = k1 + 1;
        XP1(1:2, k1) = X(1:2, i);
        MK(1, floor(i /(N + 0.1)) + 1) = MK(1, floor(i /(N+0.1)) + 1) + 1;
    end
    if ind == 2
        k2 = k2 + 1;
        XP2(1:2, k2) = X(1:2, i);
        MK(2, floor(i /(N + 0.1)) + 1) = MK(2, floor(i / (N + 0.1)) + 1) + 1;
    end
    if ind == 3
        k3 = k3 + 1;
        XP3(1:2, k3) = X(1:2, i);
        MK(3, floor(i / (N + 0.1)) + 1) = MK(3, floor(i / (N + 0.1)) + 1) + 1;
    end
end

display("Matrica konfuzije:")
display(num2str(MK))
sum(diag(MK))/sum(sum(MK))

figure,
plot(XP1(1,:), XP1(2,:), 'r*');
hold on
plot(XP2(1,:), XP2(2,:), 'b*');
plot(XP3(1,:), XP3(2,:), 'g*');
hold off