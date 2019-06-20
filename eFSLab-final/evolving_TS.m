function [centers, sigmas, conseq,R] = evolving_TS(inputs, X, radii, Omega, delta, Epsilon_min, Epsilon_max, separator, Reference_1, Reference_2, substitution_condition, creation_condition, modelo, estima, graphics)

% This function implements the algorithm "Evolving
% Takagi-Sugeno Models", proposed by P. Angelov
% 
% Inputs:   X           ... matrix with input/output data
%           modelo      ... model type
%                           0 - TS order zero
%                           1 - TS order one
%           estima      ... estimation method for consequent parameters
%                           1 - RLS (global estimation)
%                           2 - wRLS (local estimation)
%           separador   ... 3000 is the default value
% Outputs:  centers     ... centers of clusters/membership functions
%           sigmas      ... width of membership functions
%           conseq      ... parameters of rule consequents
%
% Jose Victor Ramos, August 2003
% Lara Aires, September 2008
% Matlab 6.5, Release 13
%
% Reference
% P. Angelov and D. Filev, "An Approach to On-Line Identification
% of Takagi-Sugeno Fuzzy Models"

[amostras, entradas] = size(X);
n = entradas-1;     % by default last dimension is the output
% radii = 0.3;       % usually between 0.3 and 0.5
alfa = 4/radii^2;
% Omega = 750;        % Value not too big if there is some confidence in the
                    % parameters values when covariance matrix is updated
                    % because a new rule is created
% delta = 0.1;    % grau de sobreposicao das funcoes de pertenca
                 % usually between 0.05 and 0.2

              
% =======================================================
% STAGE 1: Initialization of the rule-base structure, etc
% =======================================================

k = 1;  % from 1 to the number of samples
R = 1;  % number of fuzzy rules
Num_Regras(k) = 1;
centros(R,:) = X(k,1:n);    % cluster centers
index_rule_final(R) = k;
index_rule(R) = k;
Pontos_regras(R,:) = X(k,:);
Potential(k) = 1.0;
Potck(R) = 1.0; % Potencial of the first sample (centre)
Potential_Ref(k) = 1.0;
Potential_Min(k) = 0;
Potential_Mean(k) = 1.0;

Epsilon_down(k) = Epsilon_min*Potential_Ref(k);
Epsilon_up(k) = Epsilon_max*Potential_Ref(k);

index_rule_change = [];
index_rule_delete = [];


% Variables used in the recursive calculation of the potential for each new sample
sigmak = 0;
betak = zeros(n+1,R);
auxvk(k) = 0;
sigmakapa(k) = 0;
vek(k) = 0;


% Variables for recursive estimation of the consequent parameters
if modelo==0
    Tetak = 0;  %zeros(1,R);
    xek = 1;    %?????????????? Correcto pois neste caso as entradas so sao tidas em conta uma vez
elseif modelo ==1
    Tetak = zeros(n+1,R);       % column vector
    xek = [1, X(k,1:n)]';       % column vector xek = [ones(R,n);X(k,1:n)]';    Psik = 1*xek;   % lambda=1
end

if estima==1    % Initial conditions for RLS (estimacao global de parametros)
    if modelo==0
        Ck = Omega;
        Psik = xek;
        parametros(k,R) = Tetak';
    elseif modelo==1
        Ck = Omega*eye(R*(n+1));    % matriz quadrada (n+1)x(n+1)
        Psik = xek;  % uma vez que lambda(xk)=1, porque apenas existe uma regra
        parametros(k,:,R) = Tetak;
    end
elseif estima==2    % Initial conditions for wRLS (estimacao local de parametros)
    if modelo==0
        pik(:,:) = 0;        % column vector!! pik(R) = zeros(n+1,R);
        cik(:,:) = Omega;   
    elseif modelo==1
        pik(:,:) = zeros(n+1,R);        % column vector!! pik(R) = zeros(n+1,R);
        cik(:,:,:) = Omega*eye(n+1);    % cik(:,:,R) = Omega*eye(n+1);
    end
end

norm_cov(k) = 0;    % Atencao a maneira como esta norma e calculada!!
                    % L2-norm (maior valor singular, SVD) e nao Loo-norm (maior elemento da matriz de
                    % covariancia)

% Inicializacao da saida do modelo
saida(1) = X(1,n+1);    % since the first value that will be computed for the output
saida(2) = X(2,n+1);    % is in instant 3, i.e. output(3) WHY?

erro(1) = 0;

% Testing variables for debugging
var_teste0 = 0;
var_teste1 = 0;
var_teste2 = 0;
var_teste3 = 1; % Porque a base de regras e inicializada com uma regra
var_teste4 = 0;
var_teste5 = 0;
var_teste6 = 0;
var_teste7 = 0;
var_teste8 = 0;
var_teste9 = 0;

% Flag initialization
nova_regra = 0;

    
%%%%%%%%%%%%%%%%%%%%%%%
%   MAIN CYCLE        %
%%%%%%%%%%%%%%%%%%%%%%%

%////////////////////// CRIAÇÃO DO FICHEIRO diagnosis.txt \\\\\\\\\\\\\\\\\\\\\\\\

if exist('Diagnosis.txt', 'file') 
    delete Diagnosis.txt
end

diagnosis_file = fopen('Diagnosis.txt','a');
fprintf(diagnosis_file,'Diagnosis File \r');

%//////////////////////             FIM                  \\\\\\\\\\\\\\\\\\\\\\\\

for k=2:amostras
    
    % Update of the variables recursively computed
    % Recursive calculation of the potential of data sample (point) k
    sigmak_1 = sigmak;
    betak_1 = betak;
    
    % Update of the potential of the centers
    Potck_1 = Potck;
    
    % Update of the matrices for parameter estimation
    if estima==1    % RLS
        Ck_1 = Ck;
        Psik_1 = Psik;
    elseif estima==2    % wRLS
        pik_1 = pik;
        cik_1 = cik;
    end
    Tetak_1 = Tetak;

    % ==========================================================
    % STAGE 2: Reading the next data sample
    % ==========================================================
    
    ponto = X(k,:);     % line vector!!
    xk = ponto(1,1:n);  % inputs only!!
    yk = ponto(1,n+1);  % output
    
    % ==========================================================
    % STAGE 3: Recursive calculation of the potential of each
    %          new data sample
    %          Equation (20)
    % ==========================================================

    % Calculo do vk esquisito
    aux_vk = sum(ponto.^2);
    auxvk(k) = aux_vk;
    
    % Calculo de sigmak
    aux2 = sum(X(k-1,:).^2);
    sigmak = sigmak_1 + aux2; 
    sigmakapa(k) = sigmak;
    
    % Calculo de aux3, o v_k
    vk = 0;
    for j=1:n+1
        betak(j) = betak_1(j) + X(k-1,j);
        vk = vk + X(k,j)*betak(j);
    end
    vek(k) = vk;
    
    Potential(k) = (k-1)/((k-1)*(aux_vk+1)+sigmak-2*vk);  % Equation (20)
    % exp ??? utilizar a exponencial para calcular o potencial

    % ==========================================================
    % STAGE 4: Recursive update of the potentials of old centres
    %          taking into account the influence of the new data
    %          sample
    %          Equation (21)
    % ==========================================================

    for i=1:R
        % Computation of the distance projection between two points
        % in the axis zj (axis xj; j=1,2,...,n and axis y; j=n+1) 
        dist = 0;
        for j=1:n+1
            dist = dist + (X(k,j)-Pontos_regras(i,j)).^2;  % Actual point - cluster centre
        end
        Potck(i) = ((k-1)*Potck_1(i))/(k-2+Potck_1(i)+Potck_1(i)*dist); % Equation (21)
    end
    
    % ==========================================================
    % STAGE 5: Possible modification or upgrade of the rule base
    %          structure
    %          Equations (35), (36) and (37)
    %          Equations (27), (28) and (32)
    % ==========================================================

    Potential_Ref(k) = max(Potck);
    Potential_Min(k) = min(Potck);
    Potential_Mean(k) = sum(Potck)/R;
    Epsilon_down(k) = Epsilon_min*Potential_Ref(k);
    Epsilon_up(k) = Epsilon_max*Potential_Ref(k);
    

    distancia = 0;  % Inicializacao da variavel
    for i=1:R
        distancia(i) = sqrt(sum(abs(X(k,:)-Pontos_regras(i,:)).^2));
    end
    [Distance_Min, indice] = min(distancia);

    if Potential(k)<Potential_Min(k)
        var_teste0 = var_teste0 + 1;
    end
    if Potential(k)>Potential_Ref(k)
        var_teste1 = var_teste1 + 1;
    end

    %»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
    % CONDITION FOR SUBSTITUTION OF AN EXISTENT RULE
    %«««««««««««««««««««««««««««««««««««««««««««««««
    
    condition1=substitution_condition;
    condition2=creation_condition;    
    
    if eval(condition1)
        var_teste2 = var_teste2+1;
        nova_regra=0;
        index_rule_change(var_teste2) = k;
        
        % O novo ponto vai substituir uma regra existente
        centros(indice,:) = xk;
        index_rule_final(indice) = k;
        
        % GRANDE CHATICE!!!
        % E preciso actualizar por causa do calculo do potencial dos centros
        Pontos_regras(indice,:) = X(k,:);  
        % ISTO TEM IMPLICACOES TREMENDAS!!!
        
        Potck(indice) = Potential(k);
        
        fprintf(diagnosis_file,'\n Amostra %d Regra substituida!! %d \r', k, indice);

%         msg = sprintf('Amostra %d Regra substituida!! %d', k, indice);
%         disp(msg);
        
        if k>1 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Possivel fusao de funcoes de pertenca %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j=1:n       % coluna, ou seja entrada Xj
                for i=1:R % linha, ou seja regra Ri
                    new_c = 0;
                    new_centros_index = 0;
                    aux = 0;
                    for z=i+1:R
                        % Verifica se existem funcoes de pertenca identicas
                        if abs(centros(i,j) - centros(z,j)) < delta
                            new_c = new_c + 1;   
                            new_centros_index(new_c) = z;
                            aux = aux + centros(z,j);
                        end
                    end
                    % Actualiza os centros das funcoes resultantes da fusao
                    if new_c ~= 0
                        new_centro = (aux + centros(i,j))/(new_c+1);
                        for kapa=1:new_c
                            centros(new_centros_index(kapa),j) = new_centro;
                        end
                        centros(i,j) = new_centro;
                                                
                        var_teste8 = var_teste8+1;
                        
                        fprintf(diagnosis_file,'\n Amostra %d FUSAO (substitution rule) FP, Regra %d, Fusoes %d, Variavel X%d \r', k, i, new_c, j);
                        
%                         msg = sprintf('Amostra %d FUSAO (substitution rule) FP, Regra %d, Fusoes %d, Variavel X%d', k, i, new_c, j);
%                         disp(msg);
                    end
                end
            end
  
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Possivel fusao de regras %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % variaveis locais

            if R>1
                numR = 0;   % Para sabr quantas regras sao eliminadas
                for i=1:R
                    indices = 0;
    			    z = 0;
                    numR_local = 0;
                    indice_first = 0;
                    indices_aux = 0;
                    
                    if modelo==1
                        new_conseq = zeros(n+1,1);
                        new_psik = zeros(n+1,1);
                    elseif modelo==0
                        new_conseq = zeros(n,1);
                        new_psik = zeros(n,1);
                    end

                    new_potck = 0;
                    
                    for j=i+1:R
                        if isequal(centros(i,:),centros(j,:))   % as regras sao exactamente iguais!!
                            % eliminar regra
                            if find(indices == j)
                                % a regra ja foi adicionada
                            else
                                z = z + 1;
                                indices(z) = j;
                            end
                        end
                    end
                    indice_first = i; % Indice da regra que permanece na Base de Regras
                                    
                    % Ordenar o array de indices
					indices = sort(indices);
                    
                    % Actualiza a matriz dos centros dos clusters, a matriz de
                    % covariancia e amatriz dos parametros
                    if indices(1) ~= 0
                        numR_local = length(indices);
            			centros([indices],:)=[];
            			numR = numR + numR_local;
		
                        var_teste9 = var_teste9 + 1;
                        index_rule_delete(var_teste9) = k;  
                        
                        fprintf(diagnosis_file,'\n Amostra %d ELIMINACAO regras (substitution rule) %d \r', k, numR);
                        
%                         msg = sprintf('Amostra %d ELIMINACAO regras (substitution rule) %d', k, numR);
%                         disp(msg);
                        
                        % ACTUALIZACAO DA INFORMACAO DA REGRA i, que ficou na base de regras, as outras regras,
                        % iguais 'a i, sao eliminadas
                        indices_aux =[indice_first indices];
                        
                        % Parametros dos consequentes da regra (Piorou o comportamento!!!)
                        % Os parametros da regra que fica resultam da MEDIA dos parametros da(s) regra(s) eliminada(s)!!
                        for i=1:numR_local+1
                            if modelo==1  
                                new_conseq(:,1) = new_conseq(:,1) + Tetak((indices_aux(i)-1)*(n+1)+1:indices_aux(i)*(n+1),:);
                            elseif modelo==0
                                new_conseq(:,1) = new_conseq(:,1) + Tetak((indices_aux(i)),:);
                            end
                        end
                        
                        new_conseq(:,1) = new_conseq(:,1)/(numR_local +1);
                        
                        if modelo==0
                            Tetak(1:n,:)=new_conseq(:,1);
                        elseif modelo==1
                            Tetak((indices_aux(1)-1)*(n+1)+1:indices_aux(1)*(n+1),:) = new_conseq(:,1);
                        end
                     
                        % Matriz Psik (nao interfere, logo nao faz sentido???)
                        for i=1:numR_local+1
                            if modelo==1
                                new_psik(:,1) = new_psik(:,1) + Psik((indices_aux(i)-1)*(n+1)+1:indices_aux(i)*(n+1),:);
                            elseif modelo==0
                                new_psik(:,1) = new_psik(:,1) + Psik((indices_aux(i)),:);
                            end
                        end

                        new_psik(:,1) = new_psik(:,1)/(numR_local+1);
                        if modelo==0                         
                            Psik(1:n,:)=new_psik(:,1);
                        elseif modelo==1
                            Psik((indices_aux(1)-1)*(n+1)+1:indices_aux(1)*(n+1),:) = new_psik(:,1);
                        end
                        
                        % Potencial do centro (faz algum sentido!!!)
                        for i=1:numR_local+1
                            new_potck = new_potck + Potck(indices_aux(i));
                        end
                        new_potck = new_potck /(numR_local+1);
                        Potck(indices_aux(1)) = new_potck;
                        
                        % ELIMINACAO DA INFORMACAO DAS REGRAS REPETIDAS
                        for i=1:numR_local
                            if estima == 1
                                if modelo==1
                                    % Eliminar linhas da matriz de covariancia
                                    Ck((indices(numR_local-i+1)-1)*(n+1)+1:indices(numR_local-i+1)*(n+1),:) = [];   % n+1 linhas associadas a(s) regra(s) eliminadas
                                    % Eliminar colunas da matriz de covariancia
                                    Ck(:,(indices(numR_local-i+1)-1)*(n+1)+1:indices(numR_local-i+1)*(n+1)) = [];   % n+1 colunas associadas a(s) regra(s) eliminadas
                                
                                    % Eliminar os parametros associados a(s) regra(s) eliminadas
                                    Tetak((indices(numR_local-i+1)-1)*(n+1)+1:indices(numR_local-i+1)*(n+1),:) = [];  % Tetak = pik; mas Tetak e um vector coluna!!!
                            
                                    % Eliminar informacao associada a(s) regra(s) eliminadas
                                    Psik((indices(numR_local-i+1)-1)*(n+1)+1:indices(numR_local-i+1)*(n+1),:) = [];
                            
                                    % Eliminar a informacao associada ao potencial dos centros eliminados
                                    Potck(indices(numR_local-i+1)) = [];
                                    
                                    %Potck_1 = Potck; NAO E PRECISO!!
                                
                                elseif modelo==0                                    
                                    % Eliminar linhas da matriz de covariancia
                                    Ck((indices(numR_local-i+1)-1),:) = [];   % n+1 linhas associadas a(s) regra(s) eliminadas
                                    %Ck((indices(numR_local)),:)=[];
                                    
                                    % Eliminar colunas da matriz de covariancia
                                    Ck(:,(indices(numR_local-i+1)-1)) = [];   % n+1 colunas associadas a(s) regra(s) eliminadas
                                    %Ck(:,(indices(numR_local)))=[];                                
                                    
                                    % Eliminar os parametros associados a(s) regra(s) eliminadas
                                    Tetak((indices(numR_local-i+1)-1),:) = [];  % Tetak = pik; mas Tetak e um vector coluna!!!
                                    %Tetak((indices(numR_local)),:) = [];
                                    
                                    % Eliminar informacao associada a(s) regra(s) eliminadas
                                    Psik((indices(numR_local-i+1)-1),:) = [];
                                    %Psik((indices(numR_local)),:)=[];
                                
                                    % Eliminar a informacao associada ao potencial dos centros eliminados
                                    Potck(indices(numR_local-i+1)) = [];
                                    %Potck(indices(numR_local)) = [];
                                    
                                    %Potck_1 = Potck; NAO E PRECISO!!                                
                                end
                            end                                                                                                               
                        end
                    end
                    R = R-numR_local;                   
                    if estima == 1
                        Ck_1 = Ck;  % Para permitir a actualizacao dos parametros
                    end
                    Tetak_1 = Tetak;    % Para permitir a actualizacao dos parametros               
                
                end % for i=1:R
            end % if R>1
        end % if k>1
        
    %»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
    % CONDITION FOR CREATION OF A NEW RULE
    %«««««««««««««««««««««««««««««««««««««
    
    elseif eval(condition2)
        
        if Potential(k)<Potential_Min(k)
            var_teste5 = var_teste5 + 1;
        end
        if Potential(k)>Potential_Ref(k)
            var_teste6 = var_teste6 + 1;
        end

		var_teste3 = var_teste3+1;
		nova_regra = 1;
		R = R+1;
		centros(R,:) = xk;
        index_rule_final(R) = k;

        index_rule(var_teste3) = k;
        Pontos_regras(var_teste3,:) = X(k,:);

        % O index_rule vai dar chatices quando uma regra for eliminada!!!
        % O indice vai a vida ... pode-se mostrar no grafico as regras
        % eliminadas
        
        fprintf(diagnosis_file,'\n Amostra %d Regra criada!! %d \r', k, R);
        
%         msg = sprintf('Amostra %d Regra criada!! %d', k, R);
%         disp(msg);
        
        % E preciso alterar isto quando uma regra e eliminada!!!
        Potck(R) = Potential(k); 

        numR = 0;   % So para ter a certeza que esta tudo OK
        
        if k>1 
       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Possivel fusao de funcoes de pertenca %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j=1:n       % coluna, ou seja entrada Xj
                for i=1:R % linha, ou seja regra Ri
                    new_c = 0;
                    new_centros_index = 0;
                    aux = 0;
                    for z=i+1:R
                        % Verifica se existem funcoes de pertenca identicas
                        if abs(centros(i,j) - centros(z,j)) < delta
                            new_c = new_c + 1;   
                            new_centros_index(new_c) = z;
                            aux = aux + centros(z,j);
                        end
                    end
                    % Actualiza os centros das funcoes resultantes da fusao
                    if new_c ~= 0
                        new_centro = (aux + centros(i,j))/(new_c+1);
                        for kapa=1:new_c
                            centros(new_centros_index(kapa),j) = new_centro;
                        end
                        centros(i,j) = new_centro;
                                                
                        var_teste8 = var_teste8+1;
                        
                        fprintf(diagnosis_file,'\n Amostra %d FUSAO (new rule) FP, Regra %d, Fusoes %d, Variavel X%d \r', k, i, new_c, j);
                        
%                         msg = sprintf('Amostra %d FUSAO (new rule) FP, Regra %d, Fusoes %d, Variavel X%d', k, i, new_c, j);
%                         disp(msg);
                    end
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Possivel fusao de regras %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % A regra que foi adicionada e obrigatoriamente eliminada!!!
            
            % variaveis locais

            if R>1
                numR = 0;   % Para saber quantas regras sao eliminadas
                for i=1:R
                    indices = 0;
        			z = 0;
                    numR_local = 0;
                    
                    for j=i+1:R
                        if isequal(centros(i,:),centros(j,:))   % as regras sao exactamente iguais!!
                            % eliminar regra
                            if find(indices == j)
                                % a regra ja foi adicionada
                            else
                                z = z + 1;
                                indices(z) = j;
                            end
                        end
                    end
                
                    % Ordenar o array de indices
     			    indices = sort(indices);
                    
                    if indices(length(indices)) == R  % A regra que foi criada vai ser eliminada!!!
                  		var_teste3 = var_teste3-1;
                        index_rule(var_teste3) = [];
                        Pontos_regras(var_teste3,:) = [];
                        Potck(R) = []; 
                    end
	
     			    
                    % Vai eliminar sempre a regra mais recente, mas esta podera
                    % nao ser a melhor politica
                    
                      
                    % Actualiza a matriz dos centros dos clusters, a matriz de
                    % covariancia e amatriz dos parametros
                    if indices(1) ~= 0 %& indices(1) ~= R)
                        numR_local = length(indices);
            			centros([indices],:)=[];
            			numR = numR + numR_local;
                        
                        var_teste9 = var_teste9+1;
                        index_rule_delete(var_teste9) = k;
                        
                        fprintf(diagnosis_file,'\n\n Amostra %d ELIMINACAO regras (new rule) %d \r', k, numR);
                        
%                         msg = sprintf('Amostra %d ELIMINACAO regras (new rule) %d', k, numR);
%                         disp(msg);
                        %pause;
                        
                        % ACTUALIZACAO DA INFORMACAO DA REGRA i, que ficou
                        % na base de regras, as outras regras, iguais 'a i,
                        % sao eliminadas
                        
                        % Matriz de covariancia
                        % Parametros dos consequentes da regra                        
                        % Matriz Psik
                        % Potencial do centro
		                        
                        % ELIMINACAO DA INFORMACAO DAS REGRAS REPETIDAS
                        for i=1:numR_local-1 % numR-1, porque a informacao da regra criada nao chegou a ser adicionada                       
                            if estima==1
                                if modelo==1
                                    % Eliminar linhas da matriz de covariancia
                                    Ck((indices(numR_local-i)-1)*(n+1)+1:indices(numR_local-i)*(n+1),:) = [];   % n+1 linhas associadas a(s) regra(s) eliminadas
                                    % Eliminar colunas da matriz de covariancia
                                    Ck(:,(indices(numR_local-i)-1)*(n+1)+1:indices(numR_local-i)*(n+1)) = [];   % n+1 colunas associadas a(s) regra(s) eliminadas

                                    % Eliminar os parametros associados a(s) regra(s) eliminadas
                                    Tetak((indices(numR_local-i)-1)*(n+1)+1:indices(numR_local-i)*(n+1),:) = [];  % Tetak = pik; mas Tetak e um vector coluna!!!

                                    % Eliminar informacao associada a(s) regra(s) eliminadas
                                    Psik((indices(numR_local-i)-1)*(n+1)+1:indices(numR_local-i)*(n+1),:) = [];

                                    % Eliminar a informacao associada ao potencial dos centros eliminados
                                    Potck(indices(numR_local-i)) = [];
                                    %Potck_1 = Potck; NAO E PRECISO!!
                                elseif modelo==0
                                    % Eliminar linhas da matriz de covariancia         
                                    Ck((indices(numR_local-i)-1),:) = [];   % n+1 linhas associadas a(s) regra(s) eliminadas
                                    
                                    % Eliminar colunas da matriz de covariancia
                                    Ck(:,(indices(numR_local-i)-1)) = [];   % n+1 colunas associadas a(s) regra(s) eliminadas

                                    % Eliminar os parametros associados a(s) regra(s) eliminadas
                                    Tetak((indices(numR_local-i)-1),:) = [];  % Tetak = pik; mas Tetak e um vector coluna!!!

                                    % Eliminar informacao associada a(s) regra(s) eliminadas
                                    Psik((indices(numR_local-i)-1),:) = [];

                                    % Eliminar a informacao associada ao potencial dos centros eliminados
                                    Potck(indices(numR_local-i)) = [];                                    
                                end
                            end                          
                        end
                    end
                    R = R-numR_local;
                    
                    if estima==1
                        Ck_1 = Ck;  % Para permitir a actualizacao dos parametros
                    end
                                     
                    Tetak_1 = Tetak;    % Para permitir a actualizacao dos parametros
                 
                end % for i=1:R
            end % if R>1
        end % if k>1
        
        if numR == 0    % Se foi criada uma nova regra e nao foi eliminada de imediato

            % So no caso de ser acrescentada uma nova regra, sem ser eliminada logo de seguida, e que vai
            % executar este codigo
                        
            % «««««««««««««««««««««««««««««««««««««««««
            % Este codigo pode dar origem a uma funcao
			% »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
            
            % PORQUE E QUE ESTE CODIGO E NECESSARIO AQUI???
            % PORQUE E QUE O indice i SO VAI ATE R-1???
			sum_tau = 0;
            tau = 0;
			for i=1:R-1 % tem que ser R-1 e nao R por causa da matriz de covariancia!!
                for j=1:n
                    miu(i,j) = exp(-alfa*sum(abs(xk(1,j)-centros(i,j)).^2));  % funcoes de pertenca
                end
                tau(i) = prod(miu(i,:));    % intensidade com que dispara a regra i
                sum_tau = sum_tau + tau(i); % somatorio da intensidade com que disparam as regras
			end
            if sum_tau == 0
                
                fprintf(diagnosis_file,'\n Sum_tau e nulo, amostra =  %d \r',k);
                
%                 msg = sprintf('Sum_tau e nulo, amostra =  %d',k);
%                 disp(msg);
            end
            % Se Sum_tau e zero nenhuma regra dispara, o que nunca deveria
            % acontecer pois a saida vai ser nula!! Como se resolve???
            lambda = 0;
			for i=1:R-1 % normalizacao da intensidade com que disparam as regras
                if sum_tau == 0
                    lambda(i) = 0;
                else lambda(i) = tau(i)/sum_tau;    % a soma dos lambda's i e igual a 1
                end
			end                            
			% ««««««««««««««««««««««««««««««««««««««««««««««««««««««««««««««««««
            
            if estima==1
                if modelo==0
                    piaux = 0;
                    for i=1:R-1 % ou R?? -> tem que ser ate R-1 por causa da matriz de covariancia!!
                        piaux = piaux + lambda(i)*Tetak_1(i,1);
                    end
                    Tetak = [Tetak_1; piaux];   % acrescenta uma coluna
                    Tetak_1 = Tetak;            % Para permitir a actualizacao dos parametros
                    % ??COEFICIENTE DE ESQUECIMENTO??
                    Ro = (R^2+1)/R^2;
                    Ck = Ro*Ck_1;   % Nao tem assim tanta influencia como a partida se poderia pensar!!!
                    %Ck = Ck_1;
                    Ck(R,R) = Omega;   % n+1 linhas e colunas associadas a nova regra
                    Ck_1 = Ck;  % Para permitir a actualizacao dos parametros
                elseif modelo==1
                    piaux = zeros(n+1,1);   % vector coluna
                    for i=1:R-1 % ou R?? -> tem que ser ate R-1 por causa da matriz de covariancia!!
                        piaux = piaux + lambda(i)*Tetak_1((i-1)*(n+1)+1:i*(n+1),1);
                    end
				
                    % Reset dos parametros dos consequentes e da matriz de covariancia
                    % RLS -> Parametros dos consequentes da nova regra
                    Tetak = [Tetak_1; piaux];   % acrescenta uma coluna ???
                    Tetak_1 = Tetak;            % Para permitir a actualizacao dos parametros
				
                    % RLS -> Matriz de covariancia
                    Ro = (R^2+1)/R^2;
                    Ck = Ro*Ck_1;   % Nao tem assim tanta influencia como a partida se poderia pensar!!!
                    %Ck = Ck_1;
                    Ck((R-1)*(n+1)+1:R*(n+1),(R-1)*(n+1)+1:R*(n+1)) = eye(n+1)*Omega;   % n+1 linhas e colunas associadas a nova regra
                    Ck_1 = Ck;  % Para permitir a actualizacao dos parametros
                end
                
            elseif estima==2
                if modelo==0
                    piaux = 0;
                    for i=1:R-1
                        piaux = piaux + lambda(i)*pik_1(i);
                    end
                    pik(:,R) = piaux;   % Equacao (27a)
                    pik_1 = pik;        % Para permitir a actualizacao dos parametros
                    
                    cik(:,:)=Omega;                  
                    
                    cik_1 = cik;    % Para permitir a actualizacao dos parametros
                    for i=1:R
                        Tetak(i,1) = pik(:,i);  % Tetak = pik; mas Tetak e um vector coluna!!!
                    end
                    Tetak_1 = Tetak;    % Para permitir a actualizacao dos parametros
                elseif modelo==1
                    piaux = zeros(n+1,1);
                    for i=1:R-1
                        piaux = piaux + lambda(i)*pik_1(i);
                    end
				
                    % Reset dos parametros dos consequentes e da matriz de covariancia
                    % wRLS -> Parametros dos consequentes da nova regra
                    pik(:,R) = piaux;   % Equacao (27a)
                    pik_1 = pik;        % Para permitir a actualizacao dos parametros
                    
                    % wRLS -> Matriz de covariancia da nova regra adicionada
                    cik(:,:,R) = Omega*eye(n+1); % Equacao (32)
                    cik_1 = cik;    % Para permitir a actualizacao dos parametros
                    
                    for i=1:R
                        Tetak((i-1)*(n+1)+1:i*(n+1),1) = pik(:,i);  % Tetak = pik; mas Tetak e um vector coluna!!!
                    end
                    Tetak_1 = Tetak;    % Para permitir a actualizacao dos parametros
                end
            end
        end % if numR == 0
 
    % NO RULE IS REPLACED OR CREATED!!!
        
    else
        %nova_regra = 0;
        var_teste4 = var_teste4+1; 
    end
    
    % ===========================================================
    % STAGE 6: Recursive calculation of the consequent parameters
    %          RLS  -> Equations (24), (25) and (26)
    %          wRLS -> Equations (29), (30) and (31)
    % ===========================================================
    % This stage must be done eventhough there is a reset in parameters
    % and covariance matrices
    
    % Update of xek!!
    if modelo ==0
        %xek = 1;
        xek = [1, X(k,1:n)]';
    elseif modelo==1
        xek = [1, X(k,1:n)]';
    end
    
    % »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»
    sum_tau = 0;
    tau = 0;
    for i=1:R
        for j=1:n
            miu(i,j) = exp(-alfa*sum(abs(xk(1,j)-centros(i,j)).^2));   % membership functions
        end
        tau(i) = prod(miu(i,:));    % firing level of rule i
        sum_tau = sum_tau + tau(i); % 
    end
    if sum_tau == 0
        fprintf(diagnosis_file,'\n Sum_tau e nulo, amostra =  %d \r',k);
%         msg = sprintf('Sum_tau e nulo, amostra =  %d',k);
%         disp(msg);
    end
    % Se Sum_tau e zero nenhuma regra dispara, o que nunca deveria
    % acontecer pois a saida vai ser nula!! Como se resolve???
    lambda = 0;
    for i=1:R
        if sum_tau == 0
            lambda(i) = 0;
        else lambda(i) = tau(i)/sum_tau; % normalization of the firing level of the rules
        end
    end                             % a soma dos lambda's i e igual a 1, tens a certeza???
   % ««««««««««««««««««««««««««««««««««««««««««««««««««««««««««««««««««

    % Update of Psik!!
    if modelo ==0
        for i=1:R
            Psik(i,1) = lambda(i)*xek(1,1);
        end
    elseif modelo==1
        for i=1:R
            Psik((i-1)*(n+1)+1:i*(n+1),1) = lambda(i)*xek;
        end
    end
    
    if estima==1        % GLOBAL PARAMETER ESTIMATION (Recursive Least Squares or Kalman Filter)
        if modelo==0
            Ck = Ck_1 - (Ck_1*Psik*Psik'*Ck_1)/(1+Psik'*Ck_1*Psik); 
            Tetak = Tetak_1 + Ck*Psik*(yk-Psik'*Tetak_1);
            
        elseif modelo==1
            Ck = Ck_1 - (Ck_1*Psik*Psik'*Ck_1)/(1+Psik'*Ck_1*Psik);
            Tetak = Tetak_1 + Ck*Psik*(yk-Psik'*Tetak_1);

        end
    elseif estima==2    % LOCAL PARAMETER ESTIMATION (weighted Recursive Least Squares)
        if modelo==0
            for i=1:R             
                cik = cik_1 - ((lambda(i)*cik_1*xek'*xek*cik_1)/(1+lambda(i)*xek'*cik_1*xek)); 
            end
            for i=1:R        
                pik(i) = pik_1(i) + cik*xek'*lambda(i)*(yk-xek*pik_1(i));                          
            end
            for i=1:R
                Tetak(i,1) = pik(i);  % Tetak = pik, but Tetak is a column vector ... new parameters are added at the end!!!
            end
        elseif modelo==1
            for i=1:R
                cik(:,:,i) = cik_1(:,:,i) - (lambda(i)*cik_1(:,:,i)*xek*xek'*cik_1(:,:,i))/(1+lambda(i)*xek'*cik_1(:,:,i)*xek);
            end
            for i=1:R
                pik(:,i) = pik_1(:,i) + cik(:,:,i)*xek*lambda(i)*(yk-xek'*pik_1(:,i));
            end
            for i=1:R
                Tetak((i-1)*(n+1)+1:i*(n+1),1) = pik(:,i);  % Tetak = pik, but Tetak is a column vector ... new parameters are added at the end!!!
            end
        end
    end
    
    if estima==1
        norm_cov(k) = norm(Ck);
    elseif estima==2
        norm_cov(k) = norm(reshape(cik,[],1));
    end
    
    
    % ==========================================================
    % STAGE 7: Predicao of the output for the next time step
    %          Equation 23
    % ==========================================================

    % Prediction of model output for instant k+1
    saida(k+1) = Psik'*Tetak;   

    
    % Computation of the error between the output model and the real output
    erro(k) = saida(k)-yk;

    
    % E PRECISO ALTERAR ISTO QUANDO UMA OU MAIS REGRAS SAO ELIMINADAS
    % Tridimensional array for storing all the rules information
    if modelo==0
        for i=1:R
            parametros(k,i) = Tetak(i,1);  % Tetak is a column vector!!!     
        end        
    elseif modelo==1
        for i=1:R
            parametros(k,:,i) = Tetak((i-1)*(n+1)+1:i*(n+1),1)';  % Tetak is a column vector!!!
        end
    end
    
    fprintf(diagnosis_file,'\n Amostra %d Regras %d \r', k, R);
    
%     msg = sprintf('Amostra %d Regras %d', k, R);
%     disp(msg);

    Num_Regras(k) = R;
           
end % main cycle

fclose(diagnosis_file);

diagnosis_file = fopen('Diagnosis.txt','a');

saida= saida';

fprintf(diagnosis_file,'\r');
fprintf(diagnosis_file,'\n\n LEARNING PROCESS\r');
% msg = sprintf('\n\nLEARNING PROCESS\n');
% disp(msg);

fprintf(diagnosis_file,'\n\n Verification of condition Potential(k)<Potential_Min(k): %d \r', var_teste0);
% msg = sprintf('Verification of condition Pk(k)<Pot_inf(k): %d', var_teste0);
% disp(msg);

fprintf(diagnosis_file,'\n\n Verification of condition Potential(k)>Potential_Ref(k): %d \r', var_teste1);
% msg = sprintf('Verification of condition Pk(k)>Pot_ref(k): %d', var_teste1);
% disp(msg);

fprintf(diagnosis_file,'\n\n Number of rules created: %d\nNumber of rules replaced: %d \r', R, var_teste2);
% msg = sprintf('Number of rules created: %d\nNumber of rules replaced: %d', R, var_teste2);
% disp(msg);

fprintf(diagnosis_file,'\n\n Samples that originate new rules: %d \r');
% msg = sprintf('Samples that originate new rules: %d');
% disp(msg);

fprintf(diagnosis_file,'%d  \r', index_rule);
% msg = sprintf('%d  ', index_rule);
% disp(msg);

fprintf(diagnosis_file,'\n\n Rules created with Potential(k)<Potential_Min(k): %d \r', var_teste5);
% msg = sprintf('Rules created with Pk(k)<Pot_inf(k): %d', var_teste5);
% disp(msg);

fprintf(diagnosis_file,'\n\n Rules created with Potential(k)>Potential_Ref(k): %d \r', var_teste6);
% msg = sprintf('Rules created with Pk(k)>Pot_ref(k): %d', var_teste6);
% disp(msg);

fprintf(diagnosis_file,'\n\n Samples that originate the replacement of rules: %d \r');
% msg = sprintf('Samples that originate the replacement of rules: %d');
% disp(msg);

fprintf(diagnosis_file,'%d  \r', index_rule_change);
% msg = sprintf('%d  ', index_rule_change);
% disp(msg);


fprintf(diagnosis_file,'\n\n Samples that originate the elimination of rules: %d \r');
% msg = sprintf('Samples that originate the elimination of rules: %d');
% disp(msg);

fprintf(diagnosis_file,'%d  \r', index_rule_delete);
% msg = sprintf('%d  ', index_rule_delete);
% disp(msg);

fprintf(diagnosis_file,'\r');
fprintf(diagnosis_file,'\n\n Membership Function Simplification\r');
% msg = sprintf('\nMembership Function Simplification\n');
% disp(msg);

fprintf(diagnosis_file,'\n\n Number of modified rules that originate mf fusion: %d \r', var_teste7);
% msg = sprintf('Number of modified rules that originate mf fusion: %d', var_teste7);
% disp(msg);

fprintf(diagnosis_file,'\n\n Number of created rules that originate mf fusion: %d \r', var_teste8);
% msg = sprintf('Number of created rules that originate mf fusion: %d', var_teste8);
% disp(msg);



% Calcular o VAF (Variance Accounted For)
% Formula:
% VAF = 100%(1- var(measured data - model output)/var(measure data))
vaf=100*(1-(var(X(separator:amostras,entradas)-saida(separator:amostras))/var(X(separator:amostras,entradas))));


fprintf(diagnosis_file,'\n\n Variance Accounted For (VAF): %0.5g \r', vaf);
% msg = sprintf('\nVariance Accounted For (VAF): %0.5g', vaf);
% disp(msg);

if nargin > 2
	% Computation of performance measures for training data
	mse = (1/separator)*sum(erro(1:separator).^2);
	rmse = sqrt((1/separator)*sum(erro(1:separator).^2));
	ndei = rmse/std(X(1:separator,n+1));
    fprintf(diagnosis_file,'\n\n Performance Measures for Training:\rMSE = %0.5g \rRMSE = %0.5g \rNDEI = %0.5g \r', mse, rmse, ndei);
% 	msg = sprintf('\nPerformance Measures for Training:\n\n MSE = %0.5g \nRMSE = %0.5g \nNDEI = %0.5g \n', mse, rmse, ndei);
% 	disp(msg);
	
	% Computation of performance measures for validation data
	mse = (1/(amostras-separator))*sum(erro(separator+1:amostras).^2);
	rmse = sqrt((1/(amostras-separator))*sum(erro(separator+1:amostras).^2));
	ndei = rmse/std(X(separator+1:amostras,n+1));
    fprintf(diagnosis_file, '\n\n Performance Measures for Validation:\rMSE = %0.5g \rRMSE = %0.5g \rNDEI = %0.5g \r', mse, rmse, ndei);
% 	msg = sprintf('\nPerformance Measures for Validation:\n\n MSE = %0.5g \nRMSE = %0.5g \nNDEI = %0.5g \n', mse, rmse, ndei);
% 	disp(msg);
else
    % Calculo dos indices de desempenho para os dados de treino
	mse = (1/amostras)*sum(erro.^2);
	rmse = sqrt((1/amostras)*sum(erro.^2));
	ndei = rmse/std(X(:,n+1));
    fprintf(diagnosis_file, '\n\n Performance Measures for Training:\rMSE = %0.5g \rRMSE = %0.5g \rNDEI = %0.5g \r', mse, rmse, ndei);
% 	msg = sprintf('\nPerformance Measures for Training:\n\n MSE = %0.5g \nRMSE = %0.5g \nNDEI = %0.5g \n', mse, rmse, ndei);
% 	disp(msg);
end

% %////////////////////// FECHAR O FICHEIRO diagnosis.txt \\\\\\\\\\\\\\\\\\\\\\\\
% fclose(diagnosis_file);
% %///////////////////////////////     FIM  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

if graphics==1
    
    % Plot do grafico com a saida do modelo e do processo
    figure('Name','Model output (red) vs. real output (blue) for validation data');
    pause on;
    x = (1:amostras).';
    plot(x, saida(1:amostras), 'r', x, X(:,n+1), 'b');
    title('Model output (red) vs. real output (blue) for validation data');
    xlabel('Sample');
    ylabel('Output value');
    pause(2);
    
    % Plot do grafico com a evolucao do numero de regras
    figure('Name','Evolution of the number of rules');
    pause on;
    x = (1:amostras).';
    plot(x, Num_Regras, '-');
    title('Evolution of the number of rules');
    xlabel('Sample');
    ylabel('Rules');
    pause(2);

    % Plot do grafico com a saida do modelo
    figure('Name','Localization of rules');
    pause on;
    hold on;
    x = (1:amostras).';
    plot(x, X(:,n+1), 'b');
    title('Localization of rules (o - new rule; * - modified rule; = final rule)');
    xlabel('Sample');
    ylabel('Output value');
    
    for i=1:length(index_rule)  %R deveria coincidir com o valor de R mas ???
        value_rules(i) = X(index_rule(i),n+1);
    end
    plot(index_rule, value_rules, 'go');
    
    tam = length(index_rule_change);
    for i=1:tam
        value_rules_change(i) = X(index_rule_change(i),n+1);
    end
    if tam>0
        plot(index_rule_change, value_rules_change, 'r*');
    end
    
    tam = length(index_rule_final);
    for i=1:tam
        value_rules_final(i) = X(index_rule_final(i),n+1);
    end
    if tam>0
        plot(index_rule_final, value_rules_final, 'ks');
    end
    
    tam = length(index_rule_delete);
    for i=1:tam
        value_rules_delete(i) = X(index_rule_delete(i),n+1);
    end
    if tam >0
        plot(index_rule_delete, value_rules_delete, 'yd');
    end
    hold off;
    pause(2);
    
    % Plot do grafico com a saida real e os valores do potencial!!
    figure('Name','Potential evolution for data samples');
    x = (1:amostras).';
    %plot(x, Pk, 'g', x, Pot_ref, 'm');
    [ax, h1, h2] = plotyy(x,  X(:,n+1), x, [Potential' Potential_Min' Potential_Ref']); %mean_Potck' epson_down' epson_up']); %[Pk' Pot_ref' mean_Potck']
    title('Potential evolution for data samples');
    xlabel('Sample');
    %ylabel('Valor do potencial');
    set(get(ax(1),'Ylabel'),'String','Output value');
    set(get(ax(2),'Ylabel'),'String','Potential');
    pause(2);
    
    % Plot do grafico com a evolucao na norma da matriz de covariancia!!
    figure('Name','Norm of Covariance Matrix');
    x = (1:amostras).';
    plot(x, norm_cov);
    %plot(x, Pk, 'g', x, mean_Potck, 'm');
    title('Norm of Covariance Matrix');
    xlabel('Sample');
    ylabel('Norm value');
    pause(2);
    
    if inputs == 2
        % Plot do grafico com os clusters no espaco de entrada
        figure('Name','Clusters in input space');
        hold on;
        plot(X(:,1), X(:,2), 'o');
        for i=1:R
            plot(centros(i,1), centros(i,2), 'r*');
        end
        for i=1:R
            drawCircle(centros(i,1), centros(i,2), radii);
        end
        title('Clusters in input space');
        xlabel('First Input Values');
        ylabel('Second Input Values');
        pause(2);
        hold off;
    elseif inputs == 3
        % Plot do grafico com os clusters no espaco de entrada
        figure('Name','Clusters in input space');
        hold on;
        grid on;
        plot3(X(:,1), X(:,2), X(:,3),'o');
        for i=1:R
            plot3(centros(i,1), centros(i,2),centros(i,3), 'r*');
        end
        for i=1:R
            drawSphere(centros(i,1), centros(i,2), centros(i,3), radii);
        end
        title('Clusters in input space');
        xlabel('First Input Values');
        ylabel('Second Input Values');
        zlabel('Third Input Values');
        pause(2);
        hold off;
    else
        warndlg('Graphic representation of clusters in input space unavailable!',' Warning ');
        uiwait
    end

    
    if modelo==0
        % Plot do grafico com a evolucao dos parametros das regras
        % O numero de parametros por regra e nesta caso 1 (TS ordem 0)
        for i=1:R
            msg = sprintf('Parameters of rule %d', i);
            figure('Name',msg);
            plot(x, parametros(:,i));
            msg = sprintf('Parameters of rule %d', i);
            title(msg);
            xlabel('Samples');
            msg = sprintf('Parameters a%d0(r), a%d1(b), a%d2(g), a%d3(m), a%d4(c)', i, i, i, i, i);
            ylabel(msg);
            pause(2);
        end
        pause(2);
        hold off;
    elseif modelo==1
        % Plot do grafico com a evolucao dos parametros de uma regra
        % O numero de parametros por regra depende do numero de dimensoes
        for i=1:R
            msg = sprintf('Parameters of rule %d', i);
            figure('Name',msg);
            %plot(x, parametros(:,1,i), 'r', x, parametros(:,2,i), 'b', x, parametros(:,3,i), 'g', x, parametros(:,4,i), 'm', x, parametros(:,5,i), 'c');
            plot(x, parametros(:,1,i), 'r', x, parametros(:,2,i), 'b', x, parametros(:,3,i), 'g');
            msg = sprintf('Parameters of rule %d', i);
            title(msg);
            xlabel('Samples');
            msg = sprintf('Parameters a%d0(r), a%d1(b), a%d2(g), a%d3(m), a%d4(c)', i, i, i, i, i);
            ylabel(msg);
            pause(2);
        end
        pause(2);
        hold off;
    end
end

% Centros dos clusters
centers = centros;
fprintf(diagnosis_file, '\n\n Centers \r');
fprintf(diagnosis_file, '\n\n %d \r',centers);

% Calculo dos sigmas para os clusters
for i=1:n
    min_x(i) = min(X(:,i));
    max_x(i) = max(X(:,i));
end
sigmas = radii.*(max_x - min_x)/sqrt(8.0);

fprintf(diagnosis_file,'\n\n');
fprintf(diagnosis_file, '\n\n Sigmas \r');
fprintf(diagnosis_file, '\n\n %s \r',sigmas);

% Consequentes das regras
if modelo==0
	for i=1:R
        conseq(i) = parametros(k,i);
	end    
elseif modelo==1
	for i=1:R
        conseq(i,:) = parametros(k,:,i);
	end
end

fprintf(diagnosis_file,'\n\n');
fprintf(diagnosis_file, '\n\n Consequents \r');
fprintf(diagnosis_file, '\n\n %d \r',conseq);

%////////////////////// FECHAR O FICHEIRO diagnosis.txt \\\\\\\\\\\\\\\\\\\\\\\\
fclose(diagnosis_file);
%///////////////////////////////     FIM  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



