function tf = AveragingDFT(Avgs, umat, ymat)
    [a, b] = size(umat);
    Umat = fft(umat);
    Ymat = fft(ymat);
    U = zeros(1, a);
    Y = zeros(1, a);
    for i = 1:a
     U(i) = 1/Avgs*sum(Umat(i,:));
     Y(i) = 1/Avgs*sum(Ymat(i,:));
    end
    tf = Y./U;
end