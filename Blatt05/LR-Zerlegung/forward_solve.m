function z = forward_solve(LR,b)
	[n,n] = size(LR);
	z = zeros(n,1);
	for i=1:n
		z(i) = b(i) - sum(LR(i,1:i-1)'.*z(1:i-1)); 
	end
end
