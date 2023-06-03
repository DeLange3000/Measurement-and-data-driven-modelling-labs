function tf = AveragingAutoPowerInput(Avgs, umat, ymat)
    [a, b] = size(umat);
    Umat = fft(umat);
    Ymat = fft(ymat);
    SYU = zeros(1, a);
    SUU = zeros(1, a);
    for i = 1:a
     SYU(i) = 1/Avgs*sum(Ymat(i,:).*conj(Umat(i,:)));
     SUU(i) = 1/Avgs*sum(Umat(i,:).*conj(Umat(i,:)));
    end
    tf = SYU./SUU;
end