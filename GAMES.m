classdef GAMES < handle
    properties (SetAccess=private)
        Dim (1,1) {mustBeInteger, mustBeNonnegative}
        FuncHandle function_handle
        GridNorm (:,2) {mustBeNumeric}
        DT
        NodesNum (1,1) {mustBeInteger, mustBeNonnegative}
        ApproxNodes {mustBeNumeric}
        Tol (1,1) {mustBeLessThan(Tol, 1)}
    end

    properties (Dependent, AbortSet)
        FuncEvalLim (1,1) {mustBeInteger, mustBeGreaterThan(FuncEvalLim, 8)}
    end

    properties (SetAccess=private, Dependent)
        FuncEval (:,1) {mustBeNumeric}
    end

    properties (Access=private, Hidden)
        CandPairs
        SimplexProp
        FuncEvalStorage
        FuncEvalLimStorage
        LimNodes
    end

    methods
        function obj = GAMES(FuncHandle, RegBound, Tol, FuncEvalLim)
            arguments
                FuncHandle function_handle
                RegBound (:,2) {mustBeNumeric, mustBeReal}
                Tol (1,1) {mustBeLessThan(Tol, 1)}
                FuncEvalLim (1,1) {mustBeInteger, mustBeGreaterThan(FuncEvalLim, 8)}
            end

            obj.Dim = nargin(FuncHandle) + 1;

            if size(RegBound,1) ~= obj.Dim
                error('Invalid argument at position 2. Number of rows must match the number of function input arguments.');
            end
    
            if ~issorted(RegBound,2)
                error('Invalid argument at position 2. Rows must be presented in an ascending order.');
            end

            obj.FuncHandle = FuncHandle;
            obj.Tol = Tol;
            obj.FuncEvalLimStorage = FuncEvalLim;
            obj.FuncEvalStorage = zeros(FuncEvalLim,1);
            obj.GridNorm = zeros(obj.Dim,2);
            obj.GridNorm(:,1) = RegBound(:,1);
            obj.GridNorm(:,2) = diff(RegBound,1,2);

            if obj.Dim == 2
                [X, Y] = ndgrid([0 1]);
                InitialNodes = [X(:) Y(:)];

                obj.FuncEvalStorage(1:size(InitialNodes,1)) = obj.FuncHandle( ...
                    obj.GridNorm(1,1) + obj.GridNorm(1,2) * InitialNodes(:,1) + ...
                    1i * (obj.GridNorm(2,1) + obj.GridNorm(2,2) * InitialNodes(:,2)) ...
                    );
            else
                [X, Y, P] = ndgrid([0 1]);
                InitialNodes = [X(:) Y(:) P(:)];

                obj.FuncEvalStorage(1:size(InitialNodes,1)) = obj.FuncHandle( ...
                    obj.GridNorm(1,1) + obj.GridNorm(1,2) * InitialNodes(:,1) + ...
                    1i * (obj.GridNorm(2,1) + obj.GridNorm(2,2) * InitialNodes(:,2)), ...
                    obj.GridNorm(3,1) + obj.GridNorm(3,2) * InitialNodes(:,3) ...
                    );
            end

            obj.DT = delaunayTriangulation(InitialNodes);
            obj.NodesNum = size(obj.DT.Points,1);

            obj.SimplexProp = zeros(0,obj.Dim+2);

            obj = obj.fitTriang;
        end

        function obj = fitTriang(obj)
            if obj.Dim == 2
                obj = obj.fit2DTriang;
            else
                obj = obj.fit3DTriang;
            end
        end

        function visSol(obj, options)
            arguments
                obj
                options.VisParam (1,:) {mustBeMember(options.VisParam, [-1 1 0])} = [-1 1]
            end

            ColorSet = ...
                [
                0.8500 0.3250 0.0980;
                0 0 0;
                0 0.4470 0.7410
                ];

            if isempty(obj.ApproxNodes)
                error('No solutions were found in the given space.');
            end

            figure(Units='normalized', Position=[0.2 0.2 0.6 0.6]);

            if obj.Dim == 2
                hold on

                triplot(triangulation(obj.DT.ConnectivityList, obj.DT.Points .* obj.GridNorm(:,2).' + obj.GridNorm(:,1).'), '--k');

                for SolInd = options.VisParam
                    SolPoints = obj.ApproxNodes(obj.ApproxNodes(:,end) == SolInd,1:end-1);

                    scatter( ...
                        real(SolPoints), ...
                        imag(SolPoints), ...
                        25, ColorSet(SolInd+2,:), 'filled');
                end

                hold off
    
                xlim([obj.GridNorm(1,1) obj.GridNorm(1,1) + obj.GridNorm(1,2)]);
                ylim([obj.GridNorm(2,1) obj.GridNorm(2,1) + obj.GridNorm(2,2)]);
    
                set(gca, 'FontSize', 16);
    
                xlabel('$x$', Interpreter='latex', FontSize=18);
                ylabel('$y$', Interpreter='latex', FontSize=18);
            else
                hold on

                for SolInd = options.VisParam
                    SolPoints = obj.ApproxNodes(obj.ApproxNodes(:,end) == SolInd,1:end-1);

                    scatter3( ...
                        real(SolPoints(:,1)), ...
                        imag(SolPoints(:,1)), ...
                        SolPoints(:,2), ...
                        5, ColorSet(SolInd+2,:), 'filled');
                end

                hold off
    
                xlim([obj.GridNorm(1,1) obj.GridNorm(1,1) + obj.GridNorm(1,2)]);
                ylim([obj.GridNorm(2,1) obj.GridNorm(2,1) + obj.GridNorm(2,2)]);
                zlim([obj.GridNorm(3,1) obj.GridNorm(3,1) + obj.GridNorm(3,2)]);

                view(30,45);

                set(gca, 'FontSize', 16);
    
                xlabel('$x$', Interpreter='latex', FontSize=18);
                ylabel('$y$', Interpreter='latex', FontSize=18);
                zlabel('$p$', Interpreter='latex', FontSize=18);
            end
        end

        function value = get.FuncEval(obj)
            value = obj.FuncEvalStorage(1:obj.NodesNum);
        end

        function set.FuncEvalLim(obj, value)
            if value < obj.NodesNum
                error("Error setting property 'FuncEvalMax' of class 'triangSol'. The value must be less or equal to the number of the already existing nodes in the triangulation.");
            else
                obj.FuncEvalLimStorage = value;
            end
        end

        function value = get.FuncEvalLim(obj)
            value = obj.FuncEvalLimStorage;
        end
    end

    methods (Access=private)
        function obj = fit2DTriang(obj)
            TolSimplexNum = 0;

            if ~isempty(obj.LimNodes)
                if size(obj.LimNodes,1) > obj.FuncEvalLim - obj.NodesNum
                    NodesLim = min(size(obj.LimNodes,1), obj.FuncEvalLim - obj.NodesNum);
                    obj.LimNodes = obj.LimNodes(NodesLim+1:end,:);
                    NewNodes = obj.LimNodes(1:NodesLim,:);
                else
                    NewNodes = obj.LimNodes;
                end

                NewNodesInd = obj.NodesNum+(1:size(NewNodes,1));
                obj.DT.Points(NewNodesInd,:) = NewNodes;

                obj.FuncEvalStorage(NewNodesInd) = obj.FuncHandle( ...
                    obj.GridNorm(1,1) + obj.GridNorm(1,2) * NewNodes(:,1) + ...
                    1i * (obj.GridNorm(2,1) + obj.GridNorm(2,2) * NewNodes(:,2)) ...
                    );

                obj.NodesNum = size(obj.DT.Points,1);
                fprintf('Number of triangulation nodes: %i \n', obj.NodesNum);
            end

            while obj.NodesNum < obj.FuncEvalLim
                [CandRefNodes, CandRefSimplex, GradSimplex] = obj.getCandRefNodes;

                GradSimplex = sort(GradSimplex,2);
                [PrevInd, PrevLoc] = ismember(obj.SimplexProp(:,1:obj.Dim+1), GradSimplex, 'rows');
                PrevLoc = PrevLoc(PrevLoc ~= 0);
                obj.SimplexProp = obj.SimplexProp(PrevInd,:);
                GradSimplex(PrevLoc,:) = [];

                PrevTolSimplexNum = TolSimplexNum;

                CandRefPairs = obj.get2DSimplexPairs(CandRefSimplex);
                [GradRefPairs, TolInd] = obj.get2DSimplexPairs(obj.get2DGradRefSimplex(obj.update2DGrad(GradSimplex), TolSimplexNum));
                TolSimplexNum = sum(all(reshape(TolInd,[],2),2));

                AddRefPairs = ...
                    [
                    CandRefPairs;
                    GradRefPairs
                    ];

                AddRefNodes = (obj.DT.Points(AddRefPairs(:,1),:) + obj.DT.Points(AddRefPairs(:,2),:)) / 2;

                NewNodes = unique(...
                    [
                    CandRefNodes;
                    AddRefNodes
                    ] ...
                    , 'rows');

                if size(NewNodes,1) == 0 && PrevTolSimplexNum == TolSimplexNum
                    break
                end

                if size(NewNodes,1) > obj.FuncEvalLim - obj.NodesNum
                    NodesLim = min(size(NewNodes,1), obj.FuncEvalLim - obj.NodesNum);
                    obj.LimNodes = NewNodes(NodesLim+1:end,:);
                    NewNodes = NewNodes(1:NodesLim,:);
                end

                NewNodesInd = obj.NodesNum+(1:size(NewNodes,1));
                obj.DT.Points(NewNodesInd,:) = NewNodes;
                
                obj.FuncEvalStorage(NewNodesInd) = obj.FuncHandle( ...
                    obj.GridNorm(1,1) + obj.GridNorm(1,2) * NewNodes(:,1) + ...
                    1i * (obj.GridNorm(2,1) + obj.GridNorm(2,2) * NewNodes(:,2)) ...
                    );

                obj.NodesNum = size(obj.DT.Points,1);
                fprintf('Number of triangulation nodes: %i \n', obj.NodesNum);
            end

            Pairs = edges(obj.DT);
            obj.CandPairs = Pairs(abs(obj.argumentDiff(obj.FuncEvalStorage(Pairs(:,1)) ./ obj.FuncEvalStorage(Pairs(:,2)))) >= 2 * pi / 3,:);

            if isempty(obj.CandPairs)
                warning('No solutions were found in the given space.');
            else
                GradSimplex = obj.DT(:,:);
                CandPairNeighbors = edgeAttachments(obj.DT, obj.CandPairs);
                CandSimplexLoc = unique([CandPairNeighbors{:}]).';
                GradSimplex(CandSimplexLoc,:) = [];
    
                GradSimplex = sort(GradSimplex,2);
                [PrevInd, PrevLoc] = ismember(obj.SimplexProp(:,1:obj.Dim+1), GradSimplex, 'rows');
                PrevLoc = PrevLoc(PrevLoc ~= 0);
                obj.SimplexProp = obj.SimplexProp(PrevInd,:);
                GradSimplex(PrevLoc,:) = [];
    
                obj.update2DGrad(GradSimplex);

                CandSimplexNeighbors = neighbors(obj.DT, CandSimplexLoc);
                CandSimplexType = nan(size(CandSimplexNeighbors));

                BoundaryInd = ~(ismember(CandSimplexNeighbors, CandSimplexLoc) | isnan(CandSimplexNeighbors));
                CandSimplexLocMat = repmat(CandSimplexLoc,1,obj.Dim+1);
                [~, PropLoc] = ismember(sort(obj.DT(CandSimplexNeighbors(BoundaryInd),:),2), obj.SimplexProp(:,1:obj.Dim+1), 'rows');
                SimplexCenterLine = mean(permute(reshape(obj.DT.Points(reshape(obj.DT(CandSimplexNeighbors(BoundaryInd),:).',[],1),:).' - obj.DT.Points(reshape(obj.DT(CandSimplexLocMat(BoundaryInd),:).',[],1),:).',obj.Dim,obj.Dim+1,[]),[3 1 2]),3);
                CandSimplexType(BoundaryInd) = sign(dot(cross([SimplexCenterLine(:,1:2) ./ obj.GridNorm(1:2,2).' zeros(size(SimplexCenterLine,1),1)], [obj.SimplexProp(PropLoc,end-1:end) ./ obj.GridNorm(1:2,2).' zeros(length(PropLoc),1)],2), repmat([0 0 1],length(PropLoc),1),2));

                SimplexType = nan(length(CandSimplexLoc),1);
                [NeighborInd, NeighborLoc] = ismember(CandSimplexNeighbors, CandSimplexLoc);

                while any(isnan(SimplexType))
                    Ind = isnan(SimplexType) & ~all(isnan(CandSimplexType),2);
                    SimplexType(Ind) = sum(CandSimplexType(Ind,:),2,'omitmissing');
                    CandSimplexType(NeighborInd) = SimplexType(NeighborLoc(NeighborInd));
                end

                CandNodes = mean(permute(reshape(obj.DT.Points(reshape(obj.DT(CandSimplexLoc,:).',[],1),:).',obj.Dim,obj.Dim+1,[]),[3 1 2]),3);
                obj.ApproxNodes = ...
                    [ 
                    obj.GridNorm(1,1) + obj.GridNorm(1,2) * CandNodes(:,1) + ...
                    1i * (obj.GridNorm(2,1) + obj.GridNorm(2,2) * CandNodes(:,2)) ...
                    sign(SimplexType)
                    ];
            end
        end

        function GradLDSimplex = update2DGrad(obj, GradSimplex)
            SimplexVolMat = permute(reshape(obj.DT.Points(reshape(GradSimplex(:,2:end).',[],1),:).',obj.Dim,obj.Dim,[]),[2 1 3]) - repmat(reshape(obj.DT.Points(GradSimplex(:,1),:).',1,obj.Dim,[]),obj.Dim,1,1);
            PhaseDiffVec = permute(reshape(obj.argumentDiff(obj.FuncEvalStorage(GradSimplex(:,2:end)) ./ obj.FuncEvalStorage(GradSimplex(:,1))).',1,size(GradSimplex,2)-1,[]),[2 1 3]);

            SimplexGradXMat = SimplexVolMat;
            SimplexGradYMat = SimplexVolMat;
            SimplexGradXMat(:,1,:) = PhaseDiffVec;
            SimplexGradYMat(:,2,:) = PhaseDiffVec;

            obj.SimplexProp = ...
                [
                obj.SimplexProp;
                GradSimplex permute([obj.det2D(SimplexGradXMat) obj.det2D(SimplexGradYMat)] ./ obj.det2D(SimplexVolMat),[3 2 1])
                ];

            GradLDSimplex = ...
                [
                obj.SimplexProp(:,[1 2]);
                obj.SimplexProp(:,[1 3]);
                obj.SimplexProp(:,[2 3])
                ];
        end

        function GradRefSimplex = get2DGradRefSimplex(obj, GradLDSimplex, TolSimplexNum)
            [UniqueLDSimplex, FirstLoc] = unique(GradLDSimplex, 'rows', 'first');
            [~, LastLoc] = unique(GradLDSimplex, 'rows', 'last');
            FirstSimplexLoc = mod(FirstLoc-1,size(obj.SimplexProp,1)) + 1;
            LastSimplexLoc = mod(LastLoc-1,size(obj.SimplexProp,1)) + 1;

            LDSimplexNorm = fliplr(obj.DT.Points(UniqueLDSimplex(:,1),:) - obj.DT.Points(UniqueLDSimplex(:,2),:)) .* [-1 1];
            Flow = abs(dot((obj.SimplexProp(FirstSimplexLoc,end-1:end) - obj.SimplexProp(LastSimplexLoc,end-1:end)) ./ obj.GridNorm(1:2,2).', LDSimplexNorm ./ obj.GridNorm(1:2,2).', 2));

            [~, FlowInd] = maxk(Flow, ceil(length(unique(obj.SimplexProp(:,1:obj.Dim+1)))^(1 - 1 / obj.Dim)) + TolSimplexNum);

            GradRefSimplex = ...
                    [
                    obj.SimplexProp(FirstSimplexLoc(FlowInd),1:obj.Dim+1);
                    obj.SimplexProp(LastSimplexLoc(FlowInd),1:obj.Dim+1)
                    ];
        end

        function [Pairs, TolInd] = get2DSimplexPairs(obj, Points)
            Nodes = permute(reshape(obj.DT.Points(reshape(Points.',[],1),:).',obj.Dim,size(Points,2),[]),[3 2 1]);
            [MaxNorm, MaxInd] = max(vecnorm(circshift(Nodes,-1,2) - Nodes,2,3),[],2,'linear');
            Pairs = Points([MaxInd mod(MaxInd + size(Points,1) - 1, numel(Points)) + 1]);
            TolInd = (MaxNorm <= obj.Tol);
            Pairs(TolInd,:) = [];
        end
        
        function obj = fit3DTriang(obj)
            TolSimplexNum = 0;

            if ~isempty(obj.LimNodes)
                if size(obj.LimNodes,1) > obj.FuncEvalLim - obj.NodesNum
                    NodesLim = min(size(obj.LimNodes,1), obj.FuncEvalLim - obj.NodesNum);
                    obj.LimNodes = obj.LimNodes(NodesLim+1:end,:);
                    NewNodes = obj.LimNodes(1:NodesLim,:);
                else
                    NewNodes = obj.LimNodes;
                end

                NewNodesInd = obj.NodesNum+(1:size(NewNodes,1));
                obj.DT.Points(NewNodesInd,:) = NewNodes;

                obj.FuncEvalStorage(NewNodesInd) = obj.FuncHandle( ...
                    obj.GridNorm(1,1) + obj.GridNorm(1,2) * NewNodes(:,1) + ...
                    1i * (obj.GridNorm(2,1) + obj.GridNorm(2,2) * NewNodes(:,2)), ...
                    obj.GridNorm(3,1) + obj.GridNorm(3,2) * NewNodes(:,3) ...
                    );

                obj.NodesNum = size(obj.DT.Points,1);
                fprintf('Number of triangulation nodes: %i \n', obj.NodesNum);
            end

            while obj.NodesNum < obj.FuncEvalLim
                [CandRefNodes, CandRefSimplex, GradSimplex] = obj.getCandRefNodes;

                GradSimplex = sort(GradSimplex,2);
                [PrevInd, PrevLoc] = ismember(obj.SimplexProp(:,1:obj.Dim+1), GradSimplex, 'rows');
                PrevLoc = PrevLoc(PrevLoc ~= 0);
                obj.SimplexProp = obj.SimplexProp(PrevInd,:);
                GradSimplex(PrevLoc,:) = [];

                PrevTolSimplexNum = TolSimplexNum;

                CandRefPairs = obj.get3DSimplexPairs(CandRefSimplex);
                [GradRefPairs, TolInd] = obj.get3DSimplexPairs(obj.get3DGradRefSimplex(obj.update3DGrad(GradSimplex), TolSimplexNum));
                TolSimplexNum = sum(all(reshape(TolInd,[],2),2));

                AddRefPairs = ...
                    [
                    CandRefPairs;
                    GradRefPairs
                    ];

                NonCandRefNodes = (obj.DT.Points(AddRefPairs(:,1),:) + obj.DT.Points(AddRefPairs(:,2),:)) / 2;

                NewNodes = unique(...
                    [
                    CandRefNodes;
                    NonCandRefNodes
                    ] ...
                    , 'rows');

                if size(NewNodes,1) == 0 && PrevTolSimplexNum == TolSimplexNum
                    break
                end

                if size(NewNodes,1) > obj.FuncEvalLim - obj.NodesNum
                    NodesLim = min(size(NewNodes,1), obj.FuncEvalLim - obj.NodesNum);
                    obj.LimNodes = NewNodes(NodesLim+1:end,:);
                    NewNodes = NewNodes(1:NodesLim,:);
                end

                NewNodesInd = obj.NodesNum+(1:size(NewNodes,1));
                obj.DT.Points(NewNodesInd,:) = NewNodes;
                
                obj.FuncEvalStorage(NewNodesInd) = obj.FuncHandle( ...
                    obj.GridNorm(1,1) + obj.GridNorm(1,2) * NewNodes(:,1) + ...
                    1i * (obj.GridNorm(2,1) + obj.GridNorm(2,2) * NewNodes(:,2)), ...
                    obj.GridNorm(3,1) + obj.GridNorm(3,2) * NewNodes(:,3) ...
                    );

                obj.NodesNum = size(obj.DT.Points,1);
                fprintf('Number of triangulation nodes: %i \n', obj.NodesNum);
            end

            Pairs = edges(obj.DT);
            obj.CandPairs = Pairs(abs(obj.argumentDiff(obj.FuncEvalStorage(Pairs(:,1)) ./ obj.FuncEvalStorage(Pairs(:,2)))) >= 2 * pi / 3,:);

            if isempty(obj.CandPairs)
                warning('No solutions were found in the given space.');
            else
                GradSimplex = obj.DT(:,:);
                CandPairNeighbors = edgeAttachments(obj.DT, obj.CandPairs);
                CandSimplexLoc = unique([CandPairNeighbors{:}]).';
                GradSimplex(CandSimplexLoc,:) = [];
    
                GradSimplex = sort(GradSimplex,2);
                [PrevInd, PrevLoc] = ismember(obj.SimplexProp(:,1:obj.Dim+1), GradSimplex, 'rows');
                PrevLoc = PrevLoc(PrevLoc ~= 0);
                obj.SimplexProp = obj.SimplexProp(PrevInd,:);
                GradSimplex(PrevLoc,:) = [];
    
                obj.update3DGrad(GradSimplex);

                CandSimplexNeighbors = neighbors(obj.DT, CandSimplexLoc);
                CandSimplexType = nan(size(CandSimplexNeighbors));

                BoundaryInd = ~(ismember(CandSimplexNeighbors, CandSimplexLoc) | isnan(CandSimplexNeighbors));
                CandSimplexLocMat = repmat(CandSimplexLoc,1,obj.Dim+1);
                [~, PropLoc] = ismember(sort(obj.DT(CandSimplexNeighbors(BoundaryInd),:),2), obj.SimplexProp(:,1:obj.Dim+1), 'rows');
                SimplexCenterLine = mean(permute(reshape(obj.DT.Points(reshape(obj.DT(CandSimplexNeighbors(BoundaryInd),:).',[],1),:).' - obj.DT.Points(reshape(obj.DT(CandSimplexLocMat(BoundaryInd),:).',[],1),:).',obj.Dim,obj.Dim+1,[]),[3 1 2]),3);
                CandSimplexType(BoundaryInd) = sign(dot(cross([SimplexCenterLine(:,1:2) ./ obj.GridNorm(1:2,2).' zeros(size(SimplexCenterLine,1),1)], [obj.SimplexProp(PropLoc,end-1:end) ./ obj.GridNorm(1:2,2).' zeros(length(PropLoc),1)],2), repmat([0 0 1],length(PropLoc),1),2));

                SimplexType = nan(length(CandSimplexLoc),1);
                [NeighborInd, NeighborLoc] = ismember(CandSimplexNeighbors, CandSimplexLoc);

                while any(isnan(SimplexType))
                    Ind = isnan(SimplexType) & ~all(isnan(CandSimplexType),2);
                    SimplexType(Ind) = sum(CandSimplexType(Ind,:),2,'omitmissing');
                    CandSimplexType(NeighborInd) = SimplexType(NeighborLoc(NeighborInd));
                end

                CandNodes = mean(permute(reshape(obj.DT.Points(reshape(obj.DT(CandSimplexLoc,:).',[],1),:).',obj.Dim,obj.Dim+1,[]),[3 1 2]),3);
                obj.ApproxNodes = ...
                    [ 
                    obj.GridNorm(1,1) + obj.GridNorm(1,2) * CandNodes(:,1) + ...
                    1i * (obj.GridNorm(2,1) + obj.GridNorm(2,2) * CandNodes(:,2)) ...
                    obj.GridNorm(3,1) + obj.GridNorm(3,2) * CandNodes(:,3) ...
                    sign(SimplexType)
                    ];
            end
        end

        function GradLDSimplex = update3DGrad(obj, GradSimplex)
            SimplexVolMat = permute(reshape(obj.DT.Points(reshape(GradSimplex(:,2:end).',[],1),:).',obj.Dim,obj.Dim,[]),[2 1 3]) - repmat(reshape(obj.DT.Points(GradSimplex(:,1),:).',1,obj.Dim,[]),obj.Dim,1,1);
            PhaseDiffVec = permute(reshape(obj.argumentDiff(obj.FuncEvalStorage(GradSimplex(:,2:end)) ./ obj.FuncEvalStorage(GradSimplex(:,1))).',1,size(GradSimplex,2)-1,[]),[2 1 3]);

            SimplexGradXMat = SimplexVolMat;
            SimplexGradYMat = SimplexVolMat;
            SimplexGradXMat(:,1,:) = PhaseDiffVec;
            SimplexGradYMat(:,2,:) = PhaseDiffVec;

            obj.SimplexProp = ...
                [
                obj.SimplexProp;
                GradSimplex permute([obj.det3D(SimplexGradXMat) obj.det3D(SimplexGradYMat)] ./ obj.det3D(SimplexVolMat),[3 2 1])
                ];

            GradLDSimplex = ...
                [
                obj.SimplexProp(:,[1 2 3]);
                obj.SimplexProp(:,[1 2 4]);
                obj.SimplexProp(:,[1 3 4]);
                obj.SimplexProp(:,[2 3 4]);
                ];
        end

        function GradRefSimplex = get3DGradRefSimplex(obj, GradLDSimplex, TolSimplexNum)
            [UniqueLDSimplex, FirstLoc] = unique(GradLDSimplex, 'rows', 'first');
            [~, LastLoc] = unique(GradLDSimplex, 'rows', 'last');
            FirstSimplexLoc = mod(FirstLoc-1,size(obj.SimplexProp,1)) + 1;
            LastSimplexLoc = mod(LastLoc-1,size(obj.SimplexProp,1)) + 1;

            LDSimplexNorm = cross((obj.DT.Points(UniqueLDSimplex(:,1),:) - obj.DT.Points(UniqueLDSimplex(:,2),:)), (obj.DT.Points(UniqueLDSimplex(:,1),:) - obj.DT.Points(UniqueLDSimplex(:,3),:)),2) / 2;
            Flow = abs(dot((obj.SimplexProp(FirstSimplexLoc,end-1:end) - obj.SimplexProp(LastSimplexLoc,end-1:end)) ./ obj.GridNorm(1:2,2).', LDSimplexNorm(:,1:2) ./ obj.GridNorm(1:2,2).',2));

            [~, FlowInd] = maxk(Flow, ceil(length(unique(obj.SimplexProp(:,1:obj.Dim+1)))^(1 - 1 / obj.Dim)) + TolSimplexNum);

            GradRefSimplex = ...
                    [
                    obj.SimplexProp(FirstSimplexLoc(FlowInd),1:obj.Dim+1);
                    obj.SimplexProp(LastSimplexLoc(FlowInd),1:obj.Dim+1)
                    ];
        end

        function [Pairs, TolInd] = get3DSimplexPairs(obj, Points)
            ExtPoints = [Points Points(:,[1 3 2 4])];
            Nodes = permute(reshape(obj.DT.Points(reshape(ExtPoints.',[],1),:).',obj.Dim,size(ExtPoints,2),[]),[3 2 1]);
            [MaxNorm, MaxInd] = max(vecnorm(circshift(Nodes,-1,2) - Nodes,2,3),[],2,'linear');
            Pairs = ExtPoints([MaxInd mod(MaxInd + size(ExtPoints,1) - 1, numel(ExtPoints)) + 1]);
            TolInd = (MaxNorm <= obj.Tol);
            Pairs(TolInd,:) = [];
        end

        function [CandRefNodes, CandRefSimplex, GradSimplex] = getCandRefNodes(obj)
            Pairs = edges(obj.DT);
            obj.CandPairs = Pairs(abs(obj.argumentDiff(obj.FuncEvalStorage(Pairs(:,1)) ./ obj.FuncEvalStorage(Pairs(:,2)))) >= 2 * pi / 3,:);
            TolInd = (vecnorm(obj.DT.Points(obj.CandPairs(:,1),:) - obj.DT.Points(obj.CandPairs(:,2),:),2,2) > obj.Tol);
            CandRefNodes = (obj.DT.Points(obj.CandPairs(TolInd,1),:) + obj.DT.Points(obj.CandPairs(TolInd,2),:)) / 2;

            GradSimplex = obj.DT(:,:);
            CandSimplexLoc = edgeAttachments(obj.DT, obj.CandPairs);
            CandSimplexLoc = unique([CandSimplexLoc{:}]).';
            GradSimplex(CandSimplexLoc,:) = [];

            CandRefSimplex = obj.DT(CandSimplexLoc(sum(ismember(obj.DT(CandSimplexLoc,:), freeBoundary(obj.DT)),2) >= obj.Dim),:);
        end
    end

    methods (Access=private, Static)
        function value = det2D(M)
            value = M(1,1,:) .* M(2,2,:) - ...
                M(1,2,:) .* M(2,1,:);
        end

        function value = det3D(M)
            value = M(1,1,:) .* (M(2,2,:) .* M(3,3,:) - M(2,3,:) .* M(3,2,:)) - ...
                M(1,2,:) .* (M(2,1,:) .* M(3,3,:) - M(2,3,:) .* M(3,1,:)) + ...
                M(1,3,:) .* (M(2,1,:) .* M(3,2,:) - M(2,2,:) .* M(3,1,:));
        end

        function value = argumentDiff(arg)
            value = angle(arg);
            value(isinf(arg) | arg == 0 | isnan(arg)) = pi;
        end
    end
end

