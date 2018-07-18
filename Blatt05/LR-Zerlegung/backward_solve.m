function x = backward_solve(LR,z)
	[n,n] = size(LR);
	x = zeros(n,1);
	for i=n:-1:1
		x(i) = (z(i) - sum(LR(i,i+1:n)'.*x(i+1:n)))/LR(i,i); 
	end
end
