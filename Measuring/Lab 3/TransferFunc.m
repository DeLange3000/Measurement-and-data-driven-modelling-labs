function [H, stdH] = TransferFunc(umat, ymat, Avgs, Hfunction)
    Hfunction = str2func(Hfunction);
    [a, b] = size(umat);

    % takes first tf as true one to calculate stdH
    repetitions = floor(b/Avgs);
    if (repetitions == 1)
        disp('stdH cannot be calculated since not enough data present')
    end


    % takes first tf as true one to calculate stdH
    repetitions = floor(b/Avgs);
    if (repetitions == 1)
        disp('stdH cannot be calculated since not enough data present')
    end
    H = zeros(a, repetitions);
    stdH = zeros(a,1); 
    j = 1;
    for i = 1:repetitions
        H(:,i) = Hfunction(Avgs, umat(:,j:j+Avgs-1), ymat(:,j:j+Avgs-1));
        j = j + Avgs;
    end
    if repetitions > 1
    stdH = std(H.').';
    end
end
