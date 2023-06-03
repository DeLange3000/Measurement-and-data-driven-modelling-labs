function tf = AveragingFRF(Avgs, umat, ymat)
    [a, b] = size(umat);
    Umat = fft(umat);
    Ymat = fft(ymat);
    tf_temp = Ymat./Umat;
    tf = zeros(a, 1);
    for i = 1:a
        tf(i) = 1/Avgs*sum(tf_temp(i,:));
    end
end