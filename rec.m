function [] = rec(fig,M,ns)
X = dct2(fig);
X(abs(X)<30) = 0;
figure;
fig = idct2(X);
imagesc(fig);
colormap(gray(255));
title('Origin Image');
[row, col] = size(X);
N = row;
X_h = zeros(row, col);
Y = zeros(M, col);
for i = 1:col
    A = random('Normal',0,1,M,N);
    A = normc(A);
    n = random('Normal',0,ns,M,1);
    y = A*X(:,i) + n;
    Y(:,i) = y;
    r = y;
    x_h = zeros(N, 1);
    S_h = zeros(1,N);
    span = zeros(M, N);
    for k=1:N
        prod = A'*r;
        [~, S_h(k)] = max(abs(prod));
        span(:,k) = A(:, S_h(k));
        x_h(S_h(1:k)) = span(:, 1:k)\y;
        r = y - span(:,1:k)*x_h(S_h(1:k));
    end
    X_h(:,i) = x_h;
end
fig = idct2(Y);
figure;
imagesc(fig);
colormap(gray(255));
title('Corrupt Image');

fig_rec = idct2(X_h);
figure;
imagesc(fig_rec);
colormap(gray(255));
title(['Recovered Image','(M = ',mat2str(M),', N = ',mat2str(N),')']);
end

