function [] = NPT_0(N,loop)
prob = zeros(15, N);    % store the success probility
norm_e = zeros(15, N);  % store the average normalized error
for sm = 1:15
    for M = 1:N
        A = random('Normal',0,1,M,N);
        A = normc(A);
        cnt = 0;
        for q = 1:loop
            S = randperm(N, sm);    % store the x's non-zero entry position
            % generate the random x
            x_val = 1 + 9.*rand([sm 1]);
            x = zeros(N, 1);
            for i=1:sm
               x_val(i) = x_val(i)*(-1)^(randi(2));
            end
            for i=1:sm
               x(S(i)) = x_val(i); 
            end
            S = sort(S);    % sort it to compare more easily
            y = A*x;
            r = y;
            x_h = zeros(N, 1);  % the calculated x
            S_h = zeros(1,sm);  % the calculated non-zero entry position
            span = zeros(M, sm);    % to store chose columns of A
            % OMP
            for k=1:N
                % find the columns of A st max the inner product 
                prod = A'*r;
                [~, S_h(k)] = max(abs(prod));
                % put it into the span
                span(:,k) = A(:, S_h(k));
                % least square
                x_h(S_h(1:k)) = span(:, 1:k)\y;
                r = y - span(:,1:k)*x_h(S_h(1:k));
                e = norm(y - A*x_h);
                if (e < 1)
                   break; 
                end
            end
            norm_e(sm, M) = norm_e(sm, M) + norm(x - x_h)/norm(x);
            % count the times of success
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
fig = imagesc(prob);
colormap('gray');
colorbar;
title(['NPT (Noiseless case Success Probility, N = ',mat2str(N),')']);
xlabel('M');
ylabel('sm');
saveas(fig, ['P1SP_N',mat2str(N)]);

figure;
fig = imagesc(norm_e);
colormap('gray');
colorbar;
title(['NPT (Noiseless case Normalized Error, N = ',mat2str(N),')']);
xlabel('M');
ylabel('sm');
saveas(fig, ['P1NE_N',mat2str(N)]);
end

