classdef MOEADDAM < ALGORITHM

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [delta,nr] = Algorithm.ParameterSet(0.9,2);

            %% Generate the weight vectors
            [W,Problem.N] = UniformlyRandomlyPoint(Problem.N,Problem.M);
            % Transformation on W
            W = 1./W./repmat(sum(1./W,2),1,size(W,2));
            % Size of neighborhood
            T = ceil(Problem.N/10);
            % Size of external elite
            nEP = ceil(Problem.N*2);
            % Ratio of updated weight vectors
            nus = 0.05;
            t = 5;
            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            Z          = min(Population.objs,[],1);

            %% Optimization
            EP = [];
            adaptation_moment = round(ceil(Problem.maxFE/Problem.N)*0.05);
            tmax = 0;
            [N,M] = size(Population.objs);
            G = N;
            while Algorithm.NotTerminated(Population)
                % For each solution	
                if t > 0
                        [Population, t] = greedySelection(Population, W, Z,G)
                end
                tmax = max(t,tmax);
                T11  = 1.8 * (log(N) + sqrt(M));
                G = ceil(T11 * atan(t / tmax)  + 2);
                Offsprings(1:Problem.N) = SOLUTION();
                for i = 1 : Problem.N
                    % Choose the parents
                    if rand < delta
                        P = B(i,randperm(size(B,2)));
                    else
                        P = randperm(Problem.N);
                    end
                    % Generate an offspring
                    Offsprings(i) = OperatorGAhalf(Problem,Population(P(1:2)));
                    % Update the ideal point
                    Z = min(Z,Offsprings(i).obj);

                    g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offsprings(i).obj-Z),length(P),1).*W(P,:),[],2);
                    Population(P(find(g_old>=g_new,nr))) = Offsprings(i);
 
                end
                EP = [EP,Offsprings];
                if Problem.FE/Problem.maxFE <= 0.9
                    EP = unique(EP);
                    EP = EP(NDSort(EP.objs,1)==1);
                    if mod(ceil(Problem.FE/Problem.N),adaptation_moment) == 0
                        % Adaptive weight adjustment          
                        [Population,W] = updateWeight(Population,W,Z,EP,nus*Problem.N); 
                        B = pdist2(W,W);
                        [~,B] = sort(B,2);
                        B = B(:,1:T);
                    end
                end
            end
        end
    end
end