%% Setup
%num = 200;
opts = detectImportOptions("datos_unicos.txt");
%opts.SelectedVariableNames = [2:3];
opts.DataLines = [1, num];
M = readmatrix("datos_unicos.txt", opts);

%idks - Aristas del grafo
idxs = nchoosek(1:num, 2);
%D - Peso de cada arista
D = hypot(M(idxs(:, 1), 1) - M(idxs(:, 2), 1), ...
          M(idxs(:, 1), 2) - M(idxs(:, 2), 2));
tot_edges = length(D);

%% Constraints
% Primera restriccion: Suma de caminos activos = num de ciudades
Aeq = spones(1:length(idxs)); % Adds up the number of trips
beq = num;

% Restriccion 2: 2 aristas activas por ciudad
Aeq = [Aeq;spalloc(num,length(idxs),num*(num-1))];
% Iterar para cada ciudad
for ii = 1:num
    % Encontrar las aristas que involucran a la ciudad
    whichIdxs = (idxs == ii); 
    % Pasar a una sola columna sumando las 2 columnas
    whichIdxs = sparse(sum(whichIdxs,2)); 
    % El ii+1 renglon de A es que la suma de las aristas de A que
    % involucran a la ii-esima ciudad sea 2
    Aeq(ii+1,:) = whichIdxs';
end
beq = [beq; 2*ones(num,1)];

%% Resolver el problema inicial
intcon = 1:tot_edges;
lb = zeros(tot_edges,1);
ub = ones(tot_edges,1);

opts = optimoptions('intlinprog','Display','off');
tic;
[x_tsp,costopt,exitflag,output] = intlinprog(D,intcon,[],[],Aeq,beq,lb,ub,opts);

%num_subtours = length(detectSubtours(x_tsp, idxs))
%% Iterar hasta evitar sub recorridos
tours = detectSubtours(x_tsp,idxs);
numtours = length(tours); % number of subtours
%fprintf('# of subtours: %d\n',numtours);
iters = 1;

A = spalloc(0,tot_edges,0); % Allocate a sparse linear inequality constraint matrix
b = [];
while numtours > 1 % repeat until there is just one subtour
    % Add the subtour constraints
    b = [b;zeros(numtours,1)]; % allocate b
    A = [A;spalloc(numtours,tot_edges,num)]; % a guess at how many nonzeros to allocate
    for ii = 1:numtours
        rowIdx = size(A,1)+1; % Counter for indexing
        subTourIdx = tours{ii}; % Extract the current subtour
%         The next lines find all of the variables associated with the
%         particular subtour, then add an inequality constraint to prohibit
%         that subtour and all subtours that use those stops.
        variations = nchoosek(1:length(subTourIdx),2);
        for jj = 1:length(variations)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                       (sum(idxs==subTourIdx(variations(jj,2)),2));
            A(rowIdx,whichVar) = 1;
        end
        b(rowIdx) = length(subTourIdx)-1; % One less trip than subtour stops
    end

    % Try to optimize again
    [x_tsp,costopt,exitflag,output] = intlinprog(D,intcon,A,b,Aeq,beq,lb,ub,opts);
    
    % How many subtours this time?
    tours = detectSubtours(x_tsp,idxs);
    numtours = length(tours); % number of subtours
    iters = iters + 1;
    %fprintf('# of subtours: %d\n',numtours);
end
t = toc;