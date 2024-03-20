fileId = fopen("sample.txt", 'w');


for c = 1:C
    for d = 1:D
        for t = 1:T
            for r = 1:R
                for j = 1:nb_weeks
                    fprintf(fileId, '%d %d %d %d %f\n', c, d, t, r, j, sample(c, d, t, r, j));
                end
            end
        end
    end
end