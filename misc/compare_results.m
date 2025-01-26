function [re,cc] = tf_compare_results(L)


for ii = 1:numel(L)
    for jj = 1:numel(L)

        La = L{ii};
        Lb = L{jj};

        % e = vnorm(Lb-La,1)./vnorm(abs(La)+abs(Lb),1);
        e = vnorm(Lb-La,1)./vnorm(La+Lb,1);

        La_mc = La - mean(La);
        Lb_mc = Lb - mean(Lb);

        for kk = 1:size(Lb,2)
            c(kk) = corr(La(:,kk),Lb(:,kk)).^2;
        end

        re(ii,jj) = median(e,'omitnan');
        cc(ii,jj) = median(c,'omitnan');

    end
end

