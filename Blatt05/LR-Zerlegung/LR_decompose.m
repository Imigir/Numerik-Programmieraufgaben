function LR = LR_decompose(A)
	[n,n] = size(A);
	LR = A;
	for k=1:n-1
		LR(k+1:n,k) = LR(k+1:n,k)/LR(k,k);
		LR(k+1:n,k+1:n) = LR(k+1:n,k+1:n) - LR(k+1:n,k)*LR(k,k+1:n);
	end
end
