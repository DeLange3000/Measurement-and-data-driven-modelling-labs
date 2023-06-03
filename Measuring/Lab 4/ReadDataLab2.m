function [umat, ymat, time] = ReadDataLab2(N, Nrep, Drep, FileName)
    if(Nrep < Drep)
        disp('cannot drop more measurements then there are made')
        return
    end
    data = load(FileName);
    skip_measurements = Nrep - Drep;
    umat = zeros(N, Drep);
    ymat = zeros(N, Drep);

    for i = 1:Drep
        if(skip_measurements == 0)
            nameu = "Su";
            namey = "Sy";
        else
            nameu = "Su"+string(skip_measurements);
            namey = "Sy"+string(skip_measurements);
        end

        umat(:,i) = data.(nameu);
        ymat(:,i) = data.(namey);
        skip_measurements = skip_measurements + 1;
    end
end