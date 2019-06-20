function fismat = modelo_rsc(Xin, Xout, centers, sigmas, conseq, modelo)

% Funcao que cria um modelo difuso de Takagi-Sugeno
% de ordem 0 ou de ordem 1 a partir da informacao 
% fornecida pelo algoritmo recursivo ETS
% 
% Entradas:
%           centros ... centros dos clusters
%           sigmas  ... variancia das funcoes de pertenca
%           conseq  ... parametros dos consequentes
% Saidas:
%           modelo  ... estrutura do modelo difuso TS e respectivos parametros
%
% Jose Victor Ramos, Marco 2003
% Lara Aires, Outubro 2008
% Matlab 6.5, Release 13
%

[linhas, colunas] = size(centers);

numInp = colunas;
numOutp = 1;
numRule = linhas;

%////////////////////////////////// ALTEREI AKI \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
if modelo==0
    outEqns=[conseq(:,2:end), conseq(:,1)];
elseif modelo==1
    outEqns=[conseq(:,2:numInp+1), conseq(:,1)]';
end
%////////////////////////////////// ALTEREI AKI \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% Modifica organizacao dos parametros de cada regra, equacao y1 = k0+ k1*x1 + k2*x2 + k3*x3 (vector linha)
% para y1 = k1*x1 + k2*x2 + k3*x3 + k0 (vector coluna)

% Each column of outEqns now contains the output equation parameters
% for an output variable.  For example, if output variable y1 is given by
% the equation y1 = k1*x1 + k2*x2 + k3*x3 + k0, then column 1 of
% outEqns contains [k1 k2 k3 k0] for rule #1, followed by [k1 k2 k3 k0]
% for rule #2, etc.

%if verbose
%    disp('Creating FIS matrix...');
%end

% Find out the number of digits required for printing out the input,
% output, and rule numbers
numInDigit = floor(log10(numInp)) + 1;
numOutDigit = floor(log10(numOutp)) + 1;
numRuleDigit = floor(log10(numRule)) + 1;

% Find out the required size of the FIS matrix
numRow = 11 + (2 * (numInp + numOutp)) + (3 * (numInp + numOutp) * numRule);
numCol = numInp + numOutp + 2;      % number of columns required for the rule list  
strSize = 3 + numInDigit + numOutDigit; % size of name 'sug[numInp][numOutp]'
numCol = max(numCol,strSize);
strSize = 4 + numInDigit + numRuleDigit;    % size of 'in[numInp]mf[numRule]'
numCol = max(numCol,strSize);
strSize = 5 + numOutDigit + numRuleDigit;   % size of 'out[numOutp]fm[numRule]'
numCol = max(numCol,strSize);
numCol = max(numCol,7); % size of 'gaussmf' is 7


% Set the FIS name as 'sug[numInp][numOutp]'
theStr = sprintf('sug%g%g',numInp,numOutp);
fismat.name=theStr;
% FIS type
fismat.type = 'sugeno';
% Number of inputs and outputs
%fismat(3,1:2) = [numInp numOutp];
% Number of input membership functions
%fismat(4,1:numInp) = numRule * ones(1,numInp);
% Number of output membership functions
%fismat(5,1:numOutp) = numRule * ones(1,numOutp);
% Number of rules
%fismat(6,1) = numRule;
% Inference operators for and, or, imp, agg, and defuzz(???in this order, kliu)
fismat.andMethod = 'prod';
fismat.orMethod = 'probor';
fismat.impMethod = 'prod';
fismat.aggMethod = 'max';
fismat.defuzzMethod = 'wtaver';  % 'wtaver' ou 'wtsum'

rowIndex = 11;
% Set the input variable labels
for i=1:numInp
    theStr = sprintf('in%g',i);
    strSize = length(theStr);
    fismat.input(i).name = theStr;
end

% Set the output variable labels
for i=1:numOutp
    theStr = sprintf('out%g',i);
    strSize = length(theStr);
    fismat.output(i).name = theStr;
end

% Set the input variable ranges
%if length(xBounds) == 0
    % No data scaling range values were specified, use the actual minimum and
    % maximum values of the data.
    minX = min(Xin);
    maxX = max(Xin);
%else
%    minX = xBounds(1,1:numInp);
%    maxX = xBounds(2,1:numInp);
%end
ranges = [minX ; maxX]';
for i=1:numInp
   fismat.input(i).range = ranges(i,:);
end

% Set the output variable ranges
%if length(xBounds) == 0
    % No data scaling range values were specified, use the actual minimum and
    % maximum values of the data.
    minX = min(Xout);
    maxX = max(Xout);
%else
%    minX = xBounds(1,numInp+1:numInp+numOutp);
%    maxX = xBounds(2,numInp+1:numInp+numOutp);
%end
ranges = [minX ; maxX]';
for i=1:numOutp
   fismat.output(i).range = ranges(i,:);
end

% Set the input membership function labels
for i=1:numInp
    for j=1:numRule    
        theStr = sprintf('in%gmf%g',i,j);
        fismat.input(i).mf(j).name = theStr;
    end
end

% Set the output membership function labels
for i=1:numOutp
    for j=1:numRule       
        theStr = sprintf('out%gmf%g',i,j);
        fismat.output(i).mf(j).name = theStr;
    end
end

% Set the input membership function types
for i=1:numInp 
   for j=1:numRule
      fismat.input(i).mf(j).type = 'gaussmf';
   end   
end

% Set the output membership function types
if modelo==1
    for i=1:numOutp
        for j=1:numRule
            fismat.output(i).mf(j).type = 'linear';
        end
    end
elseif modelo==0
        for i=1:numOutp
            for j=1:numRule
                fismat.output(i).mf(j).type = 'constant';
            end
        end
end

% Set the input membership function parameters
colOfOnes = ones(numRule,1);    % a column of ones
for i=1:numInp
   for j=1:numRule
      fismat.input(i).mf(j).params = [sigmas(i) centers(j, i)];
   end
end
% Set the output membership function parameters
for i=1:numOutp
   for j=1:numRule
      %outParams = reshape(outEqns(:,i),numInp + 1,numRule);
%      fismat.output(i).mf(j).params = outParams(:,j)';
      fismat.output(i).mf(j).params = outEqns(:,j)';
        %%fismat.output(i).mf(j).params = outParams((numInp+1):,j);
   end
end

% Set the membership function pointers in the rule list
colOfEnum = [1:numRule]';

for j=1:numRule
   for i=1:numInp
      fismat.rule(j).antecedent(i)=colOfEnum(j);
   end
   for i=1:numOutp   
      fismat.rule(j).consequent(i) = colOfEnum(j);
   end
   
  % Set the antecedent operators and rule weights in the rule

   fismat.rule(j).weight=1;
   fismat.rule(j).connection=1;
end
save fismat.fis


