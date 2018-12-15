function [] = NPT_2(N,loop,ns)
prob = zeros(15, N);
norm_e = zeros(15, N);
for sm = 1:15
    for M = 1:N
        A = random('Normal',0,1,M,N);
        A = normc(A);
        cnt = 0;
        for q = 1:loop
            S = randperm(N, sm);
            x_val = 1 + 9.*rand([sm 1]);
            x = zeros(N, 1);
            for i=1:sm
               x_val(i) = x_val(i)*(-1)^(randi(2));
            end
            for i=1:sm
               x(S(i)) = x_val(i); 
            end
            S = sort(S);
            n = random('Normal',0,ns,M,1);
            n_norm = norm(n);
            y = A*x + n;
            r = y;
            x_h = zeros(N, 1);
            S_h = zeros(1,sm);
            span = zeros(M, sm);
            for k=1:N
                prod = A'*r;
                [~, S_h(k)] = max(abs(prod));
                span(:,k) = A(:, S_h(k));
                x_h(S_h(1:k)) = span(:, 1:k)\y;
                r = y - span(:,1:k)*x_h(S_h(1:k));
                e = norm(y - A*x_h);
                if (e < n_norm)
                   break; 
                end
            end
            norm_e(sm, M) = norm_e(sm, M) + norm(x - x_h)/norm(x);
            if size(S)==size(S_h)
                if (S - sort(S_h) == 0)
                   cnt = cnt + 1; 
                end
            end
        end
         prob(sm, M) = cnt/loop;
         norm_e(sm, M) = norm_e(sm, M)/loop;
    end
end
figure;
imagesc(prob);
colormap('gray');
colorbar;
title(['NPT (Noisy case Success Probility, s is not known, \sigma = ',mat2str(ns),', N = ',mat2str(N),')']);
xlabel('M');
ylabel('sm');

figure;
imagesc(norm_e);
colormap('gray');
colorbar;
title(['NPT (Noisy case Normalized Error, s is not known, \sigma = ',mat2str(ns),', N = ',mat2str(N),')']);
xlabel('M');
ylabel('sm');
end

