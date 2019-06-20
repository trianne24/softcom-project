
% figure('Name','Model output (red) vs. real output (blue) for validation data');    
% [linhas colunas]=size(output);
% 
% x = (1:linhas).';
% plot(x, output, 'r', x, inputs(:,4), 'b');
% 
% title('Model output (red) vs. real output (blue) for validation data');
% xlabel('Sample');
% ylabel('Output value');

% para duas dimensões
figure('Name','Model output (red) vs. real output (blue) for validation data');    
[linhas colunas]=size(output);

x = (1:linhas).';
plot(x, output, 'r', x, inputs(:,3), 'b');

title('Model output (red) vs. real output (blue) for validation data');
xlabel('Sample');
ylabel('Output value');