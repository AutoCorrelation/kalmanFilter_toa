classdef RSME
    properties
        
    end

    methods
        function obj = RSME()
        end

        function y = getRSME(~, estimatedPos)
            y = 0;
            for i = 1:size(estimatedPos, 3)
                for p = 1:size(estimatedPos, 2)
                    y = y + norm(estimatedPos(:, p, i) - [p; p]);
                end
            end
            y = y / (size(estimatedPos, 2) * size(estimatedPos, 3));
        end
    end


end