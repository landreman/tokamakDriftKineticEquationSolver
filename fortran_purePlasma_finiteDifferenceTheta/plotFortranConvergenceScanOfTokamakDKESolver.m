% Set the following line to whatever is appropriate for your system:
filename = 'dkeConvergence.dat';


figureOffset=1;

data = importdata(filename);

colIndex=1;

if ~strcmp(data.colheaders{colIndex}, 'Ntheta')
    error(['Column ',num2str(colIndex),' is different from what I am expecting.'])
end
colIndex = colIndex + 1;

if ~strcmp(data.colheaders{colIndex}, 'Nxi')
    error(['Column ',num2str(colIndex),' is different from what I am expecting.'])
end
colIndex = colIndex + 1;

if ~strcmp(data.colheaders{colIndex}, 'NL')
    error(['Column ',num2str(colIndex),' is different from what I am expecting.'])
end
colIndex = colIndex + 1;

if ~strcmp(data.colheaders{colIndex}, 'Nx')
    error(['Column ',num2str(colIndex),' is different from what I am expecting.'])
end
colIndex = colIndex + 1;

if ~strcmp(data.colheaders{colIndex}, 'NxPotentialsPerVth')
    error(['Column ',num2str(colIndex),' is different from what I am expecting.'])
end
colIndex = colIndex + 1;

if ~strcmp(data.colheaders{colIndex}, 'xMax')
    error(['Column ',num2str(colIndex),' is different from what I am expecting.'])
end
colIndex = colIndex + 1;

if ~strcmp(data.colheaders{colIndex}, 'q')
    error(['Column ',num2str(colIndex),' is different from what I am expecting.'])
end
colIndex = colIndex + 1;

if ~strcmp(data.colheaders{colIndex}, 'particleFlux')
    error(['Column ',num2str(colIndex),' is different from what I am expecting.'])
end
colIndex = colIndex + 1;

if ~strcmp(data.colheaders{colIndex}, 'k')
    error(['Column ',num2str(colIndex),' is different from what I am expecting.'])
end
colIndex = colIndex + 1;

if ~strcmp(data.colheaders{colIndex}, 'didItConverge')
    error(['Column ',num2str(colIndex),' is different from what I am expecting.'])
end
colIndex = colIndex + 1;

if ~strcmp(data.colheaders{colIndex}, 'elapsedTime')
    error(['Column ',num2str(colIndex),' is different from what I am expecting.'])
end
colIndex = colIndex + 1;


quantitiesToRecord = {'q','Particle flux','k','Did it converge','elapsed time'};
indicesOfQuantiesToPlot = 7:11;
linespecs = {'.-r','.-g','.-b','.-c','.-m','.-k','.-r','.-b'};
numQuantities = numel(quantitiesToRecord);

parametersToVary = {};
abscissae = {};
convergeds = {};
quantities = {};

% Check whether Ntheta was scanned:
index = 1;
[values, runIndices, scanIndices] = unique(data.data(:,index),'first');
if numel(values)>1
    parametersToVary{end+1} = 'N\theta';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data.data(:,index) == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error(['More than one value of input parameter ',num2str(index),' was repeated.'])
    end
    assert(numRepeatedValues>0)
    convergeds{end+1} = convergedValue;
    quantities{end+1} = data.data(runIndices, indicesOfQuantiesToPlot);
end

% Check whether Nxi was scanned:
index = 2;
[values, runIndices, scanIndices] = unique(data.data(:,index),'first');
if numel(values)>1
    parametersToVary{end+1} = 'N\xi';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data.data(:,index) == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error(['More than one value of input parameter ',num2str(index),' was repeated.'])
    end
    assert(numRepeatedValues>0)
    convergeds{end+1} = convergedValue;
    quantities{end+1} = data.data(runIndices, indicesOfQuantiesToPlot);
end

% Check whether NL was scanned:
index = 3;
[values, runIndices, scanIndices] = unique(data.data(:,index),'first');
if numel(values)>1
    parametersToVary{end+1} = 'NL';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data.data(:,index) == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error(['More than one value of input parameter ',num2str(index),' was repeated.'])
    end
    assert(numRepeatedValues>0)
    convergeds{end+1} = convergedValue;
    quantities{end+1} = data.data(runIndices, indicesOfQuantiesToPlot);
end

% Check whether Nx was scanned:
index = 4;
[values, runIndices, scanIndices] = unique(data.data(:,index),'first');
if numel(values)>1
    parametersToVary{end+1} = 'Nx';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data.data(:,index) == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error(['More than one value of input parameter ',num2str(index),' was repeated.'])
    end
    assert(numRepeatedValues>0)
    convergeds{end+1} = convergedValue;
    quantities{end+1} = data.data(runIndices, indicesOfQuantiesToPlot);
end

% Check whether NxPotentialsPerVth was scanned:
index = 5;
[values, runIndices, scanIndices] = unique(data.data(:,index),'first');
if numel(values)>1
    parametersToVary{end+1} = 'NxPotentialsPerVth';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data.data(:,index) == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error(['More than one value of input parameter ',num2str(index),' was repeated.'])
    end
    assert(numRepeatedValues>0)
    convergeds{end+1} = convergedValue;
    quantities{end+1} = data.data(runIndices, indicesOfQuantiesToPlot);
end

% Check whether xMax was scanned:
index = 6;
[values, runIndices, scanIndices] = unique(data.data(:,index),'first');
if numel(values)>1
    parametersToVary{end+1} = 'xMax';
    abscissae{end+1} = values;
    % Ensure only one value is repeated:
    numRepeatedValues = 0;
    convergedValue = 0;
    for i=1:numel(values)
        indices = find(data.data(:,index) == values(i));
        if numel(indices)>1
            numRepeatedValues = numRepeatedValues + 1;
            convergedValue = values(i);
        end
    end
    if numRepeatedValues>1
        error(['More than one value of input parameter ',num2str(index),' was repeated.'])
    end
    convergeds{end+1} = convergedValue;
    quantities{end+1} = data.data(runIndices, indicesOfQuantiesToPlot);
end


numParameters = numel(parametersToVary);


maxs=ones(numQuantities,1)*(-1e10);
mins=ones(numQuantities,1)*(1e10);
for iParameter = 1:numParameters
    maxs = max([maxs, quantities{iParameter}'],[],2);
    mins = min([mins, quantities{iParameter}'],[],2);
end


figure(1+figureOffset)
numRows = numQuantities;
numCols = numParameters;
clf
for iQuantity = 1:numQuantities
    if maxs(iQuantity) <= mins(iQuantity)
        maxs(iQuantity) = mins(iQuantity)+1;
    end
    for iParameter = 1:numParameters
        subplot(numRows, numCols, iParameter  + (iQuantity - 1)*numParameters)
        plot(1./abscissae{iParameter}, quantities{iParameter}(:,iQuantity)', linespecs{iQuantity})
        hold on
        plot(1./[convergeds{iParameter}, convergeds{iParameter}], [mins(iQuantity),maxs(iQuantity)],'k')
        ylim([mins(iQuantity), maxs(iQuantity)])
        xlabel(['1/',parametersToVary{iParameter}])
        ylabel(quantitiesToRecord{iQuantity})
    end
end
stringForTop=sprintf('Convergence scan from fortran local tokamak drift-kinetic equation solver');
annotation('textbox',[0 0.93 1 .07],'HorizontalAlignment','center',...
    'Interpreter','none','VerticalAlignment','bottom',...
    'FontSize',12,'LineStyle','none','String',stringForTop);

figure(8+figureOffset)
numRows = numQuantities;
numCols = numParameters;
clf
for iQuantity = 1:numQuantities
    if maxs(iQuantity) <= mins(iQuantity)
        maxs(iQuantity) = mins(iQuantity)+1;
    end
    for iParameter = 1:numParameters
        subplot(numRows, numCols, iParameter  + (iQuantity - 1)*numParameters)
        plot(abscissae{iParameter}, quantities{iParameter}(:,iQuantity)', linespecs{iQuantity})
        hold on
        plot([convergeds{iParameter}, convergeds{iParameter}], [mins(iQuantity),maxs(iQuantity)],'k')
        ylim([mins(iQuantity), maxs(iQuantity)])
        xlabel(parametersToVary{iParameter})
        ylabel(quantitiesToRecord{iQuantity})
    end
end

annotation('textbox',[0 0.93 1 .07],'HorizontalAlignment','center',...
    'Interpreter','none','VerticalAlignment','bottom',...
    'FontSize',12,'LineStyle','none','String',stringForTop);
