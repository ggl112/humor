function [Population, t] = greedySelection(Population, W, Z,G)

    % 假设Population是一个结构体数组，objs是其中的字段之一
    % 假设W是权重向量的矩阵，每一行对应一个权重向量
    % Z是参考点或理想点的向量
    t = 0;
    numWeights = size(W, 1);
    numIndividuals = length(Population);
    BestIndices = zeros(numWeights, 1); % 存储每个权重向量的最优解索引

    % 找到每个权重向量对应的最优解
    for i = 1:numWeights
        OptimalTCH = inf; % 初始化最优解的 TCH
        selectedIndices = []; % 存储已选择的个体索引
        for j = 1:numIndividuals
            currentTCH = max(abs(Population(j).objs - Z) .* W(i,:), [], 2);
            if sum(BestIndices == j) < G && currentTCH < OptimalTCH
                OptimalTCH = currentTCH;
                BestIndices(i) = j; % 更新最优解索引
                selectedIndices = [selectedIndices, j]; % 将已选择的个体索引加入列表
            end
        end
    end

    % 重组种群
    ReorganizedPopulation = Population;
    for i = 1:numWeights
        bestIndex = BestIndices(i);
        ReorganizedPopulation(i) = Population(bestIndex);
        if i ~= bestIndex
            t = t + 1;
        end
    end
    % 返回新的种群和最优解索引
    Population = ReorganizedPopulation;
end
