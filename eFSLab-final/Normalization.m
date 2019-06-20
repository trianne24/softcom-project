
[l c]= size(inputs_f);

for i=1:c
     Min(i) = min(inputs_f(:,i));
     Max(i) = max(inputs_f(:,i));
     inputs_f(:,i) = (inputs_f(:,i)-Min(i))/(Max(i) - Min(i));
 end