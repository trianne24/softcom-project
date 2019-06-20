function mam_fismat = sug2mam(fismat,sigmas,conseq,tipo_MF,part_numb_trimf, part_numb_gbellmf, slope_gbellmf, left_slope_dsigmf, part_numb_dsigmf, right_slope_dsigmf)

%SUG2MAM Transform from Sugeno to Mamdani fuzzy inference system.
% Funcao que transforma um modelo difuso de Takagi-Sugeno
% de ordem 0 ou de ordem 1 num modelo difuso de Mamdani
% 
% Entradas:
%           fismat    ... modelo difuso de Takagi-Sugeno   
%           conseq    ... parametros dos consequentes
%           tipo_MF ... tipo de função de pertença
%           n_part    ... número de partições
% Saidas:
%           modelo    ... estrutura do modelo difuso Mamdani e respectivos parametros
%
% Lara Aires, Outubro 2008
% Matlab 7.6
%

if nargin < 1,
	error('Need a FIS matrix as an input argument.');
end
if ~strcmp(fismat.type, 'sugeno'),
	error('Given FIS matrix is not a Sugeno system!');
end

in_n = length(fismat.input);
out_n = length(fismat.output);

for i = 1:out_n  
    mf_n = length(fismat.output(i).mf);
    range=fismat.output(i).range;
     x = linspace(range(1), range(2), 101);
 	for j = 1:mf_n
        fismat.output(i).mf(j).type=tipo_MF;
        
        switch tipo_MF
            case 'trimf'
                abertura=(range(2)- range(1))/(part_numb_trimf-1);
                fismat.output(i).mf(j).params=[(conseq(j)-abertura) conseq(j) (conseq(j)+abertura)];
            case 'gaussmf'
                fismat.output(i).mf(j).params=[sigmas(1) conseq(j)];
            case 'gbellmf'
                a=(range(2)-range(1))/part_numb_gbellmf;
                b=slope_gbellmf;    %b=2.5;
                c=conseq(j);
                fismat.output(i).mf(j).params=[a b c];
            case 'dsigmf'
                %assumindo que o centro da função de pertença está
                %localizado no centro correspondente dos consequentes:
                a1=left_slope_dsigmf;      %a1=63.9;
                c1=conseq(j)-((range(2)-range(1))*2/part_numb_dsigmf);
                a2=right_slope_dsigmf;     %a2=63.9;
                c2=conseq(j)+((range(2)-range(1))*2/part_numb_dsigmf);
                
                fismat.output(i).mf(j).params=[a1 c1 a2 c2];
        end
        
        mf =evalmf(x,fismat.output(i).mf(j).params,fismat.output(i).mf(j).type);
 	end
end
 
fismat.type='mamdani';
fismat.name='mam21';
fismat.defuzzMethod='centroid';
mam_fismat = fismat;
save mam_fismat.fis

