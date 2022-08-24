function res = ryser(A)

    res = 0;
    [n,n2] = size(A);

    if(n~=n2)
        error('A must be square');
    end


    % loop over all 2^n submatricies of A
    for i=1:(2^n - 1) %The last 2^n is no rows, skip that
        % use bitget to index ever possible row-subset of A
        indx = logical(bitget(i,(1:n))');

        res = res + ((-1)^(sum(indx))) * prod(sum(A(indx,:),1));
    end

    res = res * (-1)^n;

end
