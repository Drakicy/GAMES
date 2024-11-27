classdef GAMES < handle
    properties (SetAccess=private)
        %% FuncHandle: the analyzed function
        FuncHandle function_handle
        %   function handle f(v) with dim(v) >= 1

        %% Dim: dimensionality
        Dim (1,1) {mustBeInteger, mustBeNonnegative}
        %   a positive integer greater than 1

        %% GridNorm(i,j): grid normalization
        GridNorm (2,:) {mustBeNumeric}
        %   a matrix with max(i) = 2 and max(j) = Dim

        %% Tol(i): approximation tolerance level
        Tol (1,:) {mustBeLessThan(Tol, 1)}
        %   a positive vector with max(i) = Dim and values less than 1

        %% FuncEval(i): function evaluation values
        FuncEval (:,1) {mustBeNumeric}
        %   a vector with max(i) = PointNumLim

        %% DT: Delaunay triangulation
        DT
        %   a Delaunay triangulation class

        %% PointNum: triangulation points number
        PointNum (1,1) {mustBeInteger, mustBeNonnegative}
        %   a positive integer greater than 2^Dim

        %% StopCond: stop condition
        StopCond (1,1) {mustBeMember(StopCond, [0 1])} = 1
        %   either 0 or 1

        %% ApproxPoint(i,j): approximation points and their classification
        ApproxPoint (:,:) {mustBeNumeric}
        %   a matrix with max(j) = Dim
    end

    properties (Dependent, AbortSet)
        %% PointNumLim: triangulation points number limit
        PointNumLim (1,1) {mustBeInteger}
        %   a positive integer greater or equal to 2^Dim
    end

    properties (Hidden)
        %% Display: output flag
        Display (1,1) {mustBeMember(Display, [0 1])} = 1
        %   either 0 or 1
    end

    properties (Access=private, Hidden)
        %% SimplexProp(i,j): simplex properties
        SimplexProp
        %   a matrix with max(j) = Dim + 4

        PointStorage
        %   Leftover points storage

        PointNumLimStorage
        %   PointNumLim storage
    end

    methods
        function obj = GAMES(FuncHandle, RegBound, Tol, PointNumLim, options)
            %% Global Argument-based Multidimensional Equation Solver (GAMES)
            %   FuncHandle - the analyzed function, function handle f(v) 
            %                with dim(v) >= 1
            %   RegBound(i,j) - the analyzed region boundaries, a matrix 
            %                   with max(i) = dim(v) + 1 and max(j) = 2
            %   Tol(i) - approximation tolerance level, a positive vector 
            %            with max(i) = dim(v) and values less than 1
            %   PointNumLim - triangulation points number limit, a positive 
            %                 integer greater or equal to 2^(dim(v) + 1)
            %   StopCond (optional) - stop condition, either 0 or 1 (default 1)
            %   Display (optional) - output flag, either 0 or 1 (default 1)

            arguments
                FuncHandle function_handle
                RegBound (2,:) {mustBeNumeric, mustBeReal}
                Tol (1,:) {mustBePositive, mustBeLessThan(Tol, 1)}
                PointNumLim (1,1) {mustBeInteger, mustBePositive}
                options.StopCond (1,1) {mustBeMember(options.StopCond, [0 1])} = 1
                options.Display (1,1) {mustBeMember(options.Display, [0 1])} = 1
            end

            %% Setting class properties

            obj.FuncHandle = FuncHandle;

            if size(RegBound,2) < 2 || size(RegBound,2) > 3
                error('Invalid argument at position 2. Number of columns must be equal to 2 or 3.');
            end

            if RegBound(1,:) >= RegBound(2,:)
                error('Invalid argument at position 2. Column elements must be presented in an ascending order.');
            end

            obj.Dim = size(RegBound,2);

            obj.GridNorm = zeros(2,obj.Dim);
            obj.GridNorm(1,:) = RegBound(1,:);
            obj.GridNorm(2,:) = diff(RegBound,1,1);
            MinSide = min(obj.GridNorm(2,1:2));
            SideRel = obj.GridNorm(2,1:2) / MinSide;
            obj.GridNorm(2,1:2) = MinSide;

            if length(Tol) ~= obj.Dim - 1
                error('Invalid argument at position 3. Number of elements must be equal to %i.', obj.Dim - 1);
            end

            obj.Tol = [Tol(1) Tol];
            obj.Tol(1:2) = obj.Tol(1:2) ./ SideRel;

            GridVec = cell(obj.Dim,1);
            GridVec{1} = linspace(0, 1, round(SideRel(1)) + 1);
            GridVec{2} = linspace(0, 1, round(SideRel(2)) + 1);
            for i = 3:obj.Dim
                GridVec{i} = [0 1];
            end

            GridMat = cell(1,obj.Dim);
            [GridMat{:}] = ndgrid(GridVec{:});
            InitialPoint = zeros(prod(round(SideRel) + 1) * 2^(obj.Dim - 2),obj.Dim);
            for i = 1:obj.Dim
                InitialPoint(:,i) = GridMat{i}(:);
            end

            InitialPoint(:,1:2) = InitialPoint(:,1:2) .* SideRel;

            if PointNumLim <= size(InitialPoint,1)
                error('Invalid argument at position 4. Value must be greater than %i.', size(InitialPoint,1));
            end
            
            obj.PointNumLimStorage = PointNumLim;       
            obj.FuncEval(1:size(InitialPoint,1)) = obj.FuncHandle( ...
                [
                obj.GridNorm(1,1) + obj.GridNorm(2,1) * InitialPoint(:,1) + ...
                1i * (obj.GridNorm(1,2) + obj.GridNorm(2,2) * InitialPoint(:,2)) ...
                obj.GridNorm(1,3:end) + obj.GridNorm(2,3:end) .* InitialPoint(:,3:end)
                ] ...
                );

            obj.DT = delaunayTriangulation(InitialPoint);
            obj.PointNum = size(InitialPoint,1);

            obj.SimplexProp = zeros(0,obj.Dim+4);

            obj.StopCond = options.StopCond;
            obj.Display = options.Display;

            %% Initializing triangulation fitting

            obj.fitTriang;
        end

        function obj = fitTriang(obj)
            %% fitTriang: triangulation fitting

            if ~isempty(obj.PointStorage)
                %% Completing previous run
                PointLim = min(size(obj.PointStorage,1), obj.PointNumLim - obj.PointNum);
                NewPoint = obj.PointStorage(1:PointLim,:);
                obj.PointStorage = obj.PointStorage(PointLim+1:size(obj.PointStorage,1),:);

                obj.updateDT(NewPoint);

                %% Displaying algorithm progress

                if obj.Display
                    fprintf('Number of triangulation points: %i\n', obj.PointNum);
                end
            end

            while obj.PointNum < obj.PointNumLim
                %% Identifying candidate refinement pairs

                Pair = edges(obj.DT);
                CandPair = Pair(abs(obj.argDiff(Pair(:,1), Pair(:,2))) >= 2 * pi / 3,:);
                CandSimplexLoc = edgeAttachments(obj.DT, CandPair);
                CandSimplexLoc = unique([CandSimplexLoc{:}])';

                %% Updating absolute value gradient
                
                obj.updateSimplexProp(CandSimplexLoc);

                %% Identifing gradient simplex pairs

                LDSimplex = zeros(size(obj.SimplexProp,1)*(obj.Dim+1),obj.Dim);
                LDSimplexInd = nchoosek(1:obj.Dim+1,obj.Dim);
                for i = 1:obj.Dim+1
                    LDSimplex((1:size(obj.SimplexProp,1))+(i-1)*size(obj.SimplexProp,1),:) = obj.SimplexProp(:,LDSimplexInd(i,:));
                end

                [~, FirstLDSimplexLoc, UniqueLDSimplexLoc] = unique(LDSimplex, 'rows');
                SecondLDSimplexLoc = setdiff(1:size(LDSimplex,1),FirstLDSimplexLoc);
                FirstLDSimplexLoc = FirstLDSimplexLoc(UniqueLDSimplexLoc(SecondLDSimplexLoc));
                UniqueLDSimplex = LDSimplex(SecondLDSimplexLoc,:);
                FirstSimplexLoc = mod(FirstLDSimplexLoc-1,size(obj.SimplexProp,1)) + 1;
                SecondSimplexLoc = mod(SecondLDSimplexLoc-1,size(obj.SimplexProp,1)) + 1;

                %% Verifying simplex tolerance condition

                TolInd = ~(obj.SimplexProp(FirstSimplexLoc,end-2) & obj.SimplexProp(SecondSimplexLoc,end-2) & ~obj.StopCond);
                UniqueLDSimplex = UniqueLDSimplex(TolInd,:);
                FirstSimplexLoc = FirstSimplexLoc(TolInd);
                SecondSimplexLoc = SecondSimplexLoc(TolInd);

                %% Ranking gradient simplex pairs

                LDSimplexVec = permute(reshape(obj.DT.Points(reshape(UniqueLDSimplex(:,2:end)',[],1),:)',obj.Dim,obj.Dim-1,[]),[2 1 3]) - ...
                    permute(obj.DT.Points(UniqueLDSimplex(:,1),:),[3 2 1]);
                AbsFlow =  ...
                    abs( ...
                    (obj.SimplexProp(FirstSimplexLoc,end-1) - obj.SimplexProp(SecondSimplexLoc,end-1)) .* obj.detND(LDSimplexVec(:,1:end ~= 1,:)) ./  obj.GridNorm(2,1) - ...
                    (obj.SimplexProp(FirstSimplexLoc,end) - obj.SimplexProp(SecondSimplexLoc,end)) .* obj.detND(LDSimplexVec(:,1:end ~= 2,:)) ./  obj.GridNorm(2,2) ...
                    ) * prod(obj.GridNorm(2,:));

                % Absolute value refinement condition
                RefInd = (AbsFlow >= max(AbsFlow) / 2);
                % can be changed according to a certain method
                
                GradRefSimplex = unique( ...
                    [
                    obj.SimplexProp(FirstSimplexLoc(RefInd),1:obj.Dim+1);
                    obj.SimplexProp(SecondSimplexLoc(RefInd),1:obj.Dim+1)
                    ], ...
                    'rows');

                %% Identifying gradient refinement pairs

                PermuteInd = nchoosek(1:obj.Dim+1,2);
                
                GradPair = GradRefSimplex(:,PermuteInd(1,:));
                GradPairProp = abs(obj.absDiff(GradRefSimplex(:,PermuteInd(1,1)), GradRefSimplex(:,PermuteInd(1,2))));  
                for i = 2:obj.Dim+1
                    temp = abs(obj.absDiff(GradRefSimplex(:,PermuteInd(i,1)), GradRefSimplex(:,PermuteInd(i,2))));
                    PropInd = (temp > GradPairProp);
                    GradPair(PropInd,:) = GradRefSimplex(PropInd,PermuteInd(i,:));
                    GradPairProp(PropInd) = temp(PropInd);
                end

                %% Verifyng nonboundary pair tolerance condition

                CandTolInd = any(abs(obj.DT.Points(CandPair(:,1),:) - obj.DT.Points(CandPair(:,2),:)) > obj.Tol,2);
                CandPair = CandPair(CandTolInd,:);

                GradTolInd = any(abs(obj.DT.Points(GradPair(:,1),:) - obj.DT.Points(GradPair(:,2),:)) > obj.Tol,2);
                obj.SimplexProp(ismember(obj.SimplexProp(:,1:obj.Dim+1), GradRefSimplex(~GradTolInd,:), 'rows'),end-2) = 1;
                GradPair = GradPair(GradTolInd,:);

                %% Identifying boundary refinement pairs

                BoundRefSimplexLoc = edgeAttachments(obj.DT, ...
                    [
                    CandPair;
                    GradPair
                    ] ...
                    );
                BoundRefSimplexLoc = [BoundRefSimplexLoc{:}]';
                BoundRefSimplex = sort(obj.DT(BoundRefSimplexLoc,:),2);
                BoundRefSimplex = BoundRefSimplex(sum(ismember(BoundRefSimplex, freeBoundary(obj.DT)),2) == obj.Dim,:);

                BoundPair = BoundRefSimplex(:,PermuteInd(1,:));
                BoundPairDist = vecnorm(diff(permute(reshape(obj.DT.Points(reshape(BoundRefSimplex(:,PermuteInd(1,:))',[],1),:)',obj.Dim,2,[]),[3 1 2]),1,3),2,2);
                for i = 2:obj.Dim+1
                    temp = vecnorm(diff(permute(reshape(obj.DT.Points(reshape(BoundRefSimplex(:,PermuteInd(i,:))',[],1),:)',obj.Dim,2,[]),[3 1 2]),1,3),2,2);
                    DistInd = (temp > BoundPairDist);
                    BoundPair(DistInd,:) = BoundRefSimplex(DistInd,PermuteInd(i,:));
                    BoundPairDist(DistInd) = temp(DistInd);
                end

                %% Verifyng boundary pair tolerance condition

                BoundTolInd = any(abs(obj.DT.Points(BoundPair(:,1),:) - obj.DT.Points(BoundPair(:,2),:)) > obj.Tol,2);
                BoundPair = BoundPair(BoundTolInd,:);

                %% Adding new triangulation points

                NewPair = ...
                    [
                    CandPair;
                    GradPair;
                    BoundPair
                    ];
                NewPoint = unique((obj.DT.Points(NewPair(:,1),:) + obj.DT.Points(NewPair(:,2),:)) / 2, 'rows');

                if isempty(NewPoint) && obj.StopCond
                    break
                end

                PointLim = min(size(NewPoint,1), obj.PointNumLim - obj.PointNum);
                obj.PointStorage = NewPoint(PointLim+1:size(NewPoint,1),:);
                NewPoint = NewPoint(1:PointLim,:);

                obj.updateDT(NewPoint);

                %% Displaying algorithm progress

                if obj.Display
                    fprintf('Number of triangulation points: %i\n', obj.PointNum);
                end
            end

            %% Identifying candidate pairs

            Pair = edges(obj.DT);
            CandPair = Pair(abs(obj.argDiff(Pair(:,1), Pair(:,2))) >= 2 * pi / 3,:);

            %% Identifying and classifying approximation points

            if isempty(CandPair)
                warning('No solutions were found in the domain.');
            else
                %% Identifying candidate region

                CandSimplexLoc = edgeAttachments(obj.DT, CandPair);
                CandSimplexLoc = unique([CandSimplexLoc{:}])';

                %% Updating absolute value gradient

                obj.updateSimplexProp(CandSimplexLoc);

                %% Identifying boundary of the candidate region

                CandSimplexNeighbors = neighbors(obj.DT, CandSimplexLoc);
                CandSimplexType = nan(size(CandSimplexNeighbors));
                [NeighborInd, NeighborLoc] = ismember(CandSimplexNeighbors, CandSimplexLoc);
                
                %% Classifying noncandidate simplexes

                BoundInd = ~(NeighborInd | isnan(CandSimplexNeighbors));
                CandSimplexLocMat = repmat(CandSimplexLoc,1,obj.Dim+1);
                [~, PropLoc] = ismember(sort(obj.DT(CandSimplexNeighbors(BoundInd),:),2), obj.SimplexProp(:,1:obj.Dim+1), 'rows');
                CandSimplexType(BoundInd) = ...
                    dot( ...
                    mean(permute(reshape(obj.DT.Points(reshape(obj.DT(CandSimplexNeighbors(BoundInd),:)',[],1),1:2)' - obj.DT.Points(reshape(obj.DT(CandSimplexLocMat(BoundInd),:)',[],1),1:2)',2,obj.Dim+1,[]),[3 1 2]),3) .* obj.GridNorm(2,1:2), ...
                    obj.SimplexProp(PropLoc,end-1:end), ...
                    2);

                %% Classifying candidate simplexes

                SimplexType = nan(length(CandSimplexLoc),1);

                while any(isnan(SimplexType))
                    Ind = isnan(SimplexType) & ~all(isnan(CandSimplexType),2);
                    SimplexType(Ind) = mean(CandSimplexType(Ind,:),2,'omitmissing');
                    CandSimplexType(NeighborInd) = SimplexType(NeighborLoc(NeighborInd));
                end

                %% Finalizing approximation process

                TriangPoint = obj.GridNorm(1,:) + obj.GridNorm(2,:) .* mean(permute(reshape(obj.DT.Points(reshape(obj.DT(CandSimplexLoc,:)',[],1),:)',obj.Dim,obj.Dim+1,[]),[3 1 2]),3);
                obj.ApproxPoint = ...
                    [ 
                    TriangPoint(:,1) + ...
                    1i * TriangPoint(:,2) ...
                    TriangPoint(:,3:end) ...
                    sign(SimplexType)
                    ];
            end
        end

        function visTriang(obj, options)
            %% visTriang: triangulation visualization
            %   PointClass (optional) - approximation point class, a member
            %                          of [0, -1, 1]

            arguments
                obj
                options.PointClass (1,:) {mustBeMember(options.PointClass, [0 -1 1])} = [-1 1]
            end

            %% Checking for solutions

            if isempty(obj.ApproxPoint)
                warning('No solutions were found in the domain.');
                return
            end

            %% Setting figure parameters

            ColorSet = ...
                [
                0.8500 0.3250 0.0980;
                0 0 0;
                0 0.4470 0.7410
                ];

            figure(Units='normalized', Position=[0.2 0.2 0.6 0.6]);

            %% Plotting triangulation

            if obj.Dim == 2
                hold on

                TR = triangulation(obj.DT(:,:), obj.GridNorm(1,:) + obj.GridNorm(2,:) .* obj.DT.Points);
                triplot(TR, '-k');

                for SolInd = options.PointClass
                    SolPoint = obj.ApproxPoint(obj.ApproxPoint(:,end) == SolInd,1:obj.Dim-1);

                    scatter( ...
                        real(SolPoint), ...
                        imag(SolPoint), ...
                        25, ColorSet(SolInd+2,:), 'filled');
                end

                hold off
            else
                hold on

                for SolInd = options.PointClass
                    SolPoint = obj.ApproxPoint(obj.ApproxPoint(:,end) == SolInd,1:obj.Dim-1);

                    scatter3( ...
                        real(SolPoint(:,1)), ...
                        imag(SolPoint(:,1)), ...
                        SolPoint(:,2), ...
                        5, ColorSet(SolInd+2,:), 'filled');
                end

                hold off

                view(30,45);

                zlabel('$p$', Interpreter='latex', FontSize=18);
            end

            LowerLim = min(obj.GridNorm(1,:) + obj.GridNorm(2,:) .* obj.DT.Points,[],1);
            UpperLim = max(obj.GridNorm(1,:) + obj.GridNorm(2,:) .* obj.DT.Points,[],1);
            axis(reshape([LowerLim; UpperLim],[],1)');

            set(gca, 'FontSize', 16);
    
            xlabel('$x$', Interpreter='latex', FontSize=18);
            ylabel('$y$', Interpreter='latex', FontSize=18);
        end

        function set.PointNumLim(obj, value)
            if value <= obj.PointNum
                error("Error setting property 'PointNumLim' of class 'GAMES'. The value must be greater than %i.", obj.PointNum);
            else
                obj.PointNumLimStorage = value;
            end
        end

        function value = get.PointNumLim(obj)
            value = obj.PointNumLimStorage;
        end
    end

    methods (Access=private)
        function updateDT(obj, NewPoint)
            %% updateDT: Delaunay triangulation update
            %   NewPoint(i,j) - an array with max(j) = Dim

            obj.DT.Points(end+1:end+size(NewPoint,1),:) = NewPoint;
            obj.FuncEval(end+1:end+size(NewPoint,1)) = obj.FuncHandle( ...
                [
                obj.GridNorm(1,1) + obj.GridNorm(2,1) * NewPoint(:,1) + ...
                1i * (obj.GridNorm(1,2) + obj.GridNorm(2,2) * NewPoint(:,2)) ...
                obj.GridNorm(1,3:end) + obj.GridNorm(2,3:end) .* NewPoint(:,3:end)
                ] ...
                );

            obj.PointNum = size(obj.DT.Points,1);
        end

        function updateSimplexProp(obj, NonGradSimplexLoc)
            %% updateSimplexProp: Simplex properties update
            %   CandSimplexLoc(i) - a positive integer vector

            %% Constructing noncandidate region

            GradSimplex = obj.DT(:,:);
            GradSimplex(NonGradSimplexLoc,:) = [];
            GradSimplex = sort(GradSimplex,2);
            [PrevInd, PrevLoc] = ismember(obj.SimplexProp(:,1:obj.Dim+1), GradSimplex, 'rows');
            obj.SimplexProp = obj.SimplexProp(PrevInd,:);
            GradSimplex(PrevLoc(PrevInd),:) = [];

            %% Calculating absolute value gradient

            SimplexVec = permute(reshape(obj.DT.Points(reshape(GradSimplex(:,2:end)',[],1),:)',obj.Dim,obj.Dim,[]),[2 1 3]) - ...
                permute(obj.DT.Points(GradSimplex(:,1),:),[3 2 1]);
            AbsDiffVec = permute(obj.absDiff(GradSimplex(:,1), GradSimplex(:,2:end)),[2 3 1]);

            AbsVecX = SimplexVec;
            AbsVecY = SimplexVec;
            AbsVecX(:,1,:) = AbsDiffVec;
            AbsVecY(:,2,:) = AbsDiffVec;

            obj.SimplexProp = ...
                [
                obj.SimplexProp;
                GradSimplex zeros(size(GradSimplex,1),1) [obj.detND(AbsVecX) obj.detND(AbsVecY)] ./ obj.GridNorm(2,[1 2]) ./ obj.detND(SimplexVec)
                ];
        end

        function value = absDiff(obj, FirstPoint, SecondPoint)
            %% absDiff: absolute value difference between a pair of points
            %   FirstPoint(i) - a positive integer vector
            %   SecondPoint(i) - a positive integer vector

            value = log(abs(obj.FuncEval(SecondPoint) ./ obj.FuncEval(FirstPoint)));
        end

        function value = argDiff(obj, FirstPoint, SecondPoint)
            %% argDiff: argument difference between a pair of points
            %   FirstPoint(i) - a positive integer vector
            %   SecondPoint(i) - a positive integer vector

            arg = obj.FuncEval(SecondPoint) ./ obj.FuncEval(FirstPoint);
            value = angle(arg);
            value(isinf(arg) | arg == 0 | isnan(arg)) = pi;
        end

        function value = detND(obj, M, varargin)
            %% detND: vectorized matrix determinant
            %   M(i,j,k) - an array with max(i) = max(j)

            if prod(size(M,[1 2])) == 1
                value = M;
            else
                value = 0;
                for i = 1:size(M,1)
                    value = value + (-1)^(i+1) * M(i,1,:) .* obj.detND(M(1:size(M,1) ~= i,2:end,:),0);
                end
            end

            if isempty(varargin)
                value = permute(value,[3 1 2]);
            end
        end
    end
end

