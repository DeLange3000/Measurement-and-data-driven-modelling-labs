function tf = AveragingAutoPowerOutput(Avgs, umat, ymat)
    [a, b] = size(umat);
    Umat = fft(umat);
    Ymat = fft(ymat);
    SYY = zeros(1, a);
    SUY = zeros(1, a);
    for i = 1:a
     SYY(i) = 1/Avgs*sum(Ymat(i,:).*conj(Ymat(i,:)));
     SUY(i) = 1/Avgs*sum(Umat(i,:).*conj(Ymat(i,:)));
    end
    tf = SYY./SUY;
end