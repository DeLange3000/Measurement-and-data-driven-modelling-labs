function tf = AvgTimeDomain(Avgs, umat, ymat)
    [a, b] = size(umat);
    u = zeros(1, a);
    y = zeros(1, a);
    for i = 1:a
     u(i) = 1/Avgs*sum(umat(i,:));
     y(i) = 1/Avgs*sum(ymat(i,:));
    end
    tf = fft(y)./fft(u);
end