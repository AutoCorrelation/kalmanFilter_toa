function scatterParticle(input)
for k=1:size(input,4)-1
    subplot(2,2,k)
    for i=1:size(input,2)
        % for j=1:size(input,3)
            scatter(squeeze(input(1,i,:,k)),squeeze(input(2,i,:,k)));
        % end
    end
end
end