using LinearAlgebra, Distributions, HypothesisTests

function McCullagh_test(mat::Matrix{Int})
	r, c = size(mat)
    r == c || throw(DimensionMismatch("input matrix 'mat' should be a square matrix."))
    #Symmetric matrix N, the first equation in p.450
    N = zeros(Int64, r -1, c -1)
	for i = 1:r-1
		for j = i:r-1
			N[i,j] = N[j, i] = sum(mat[1:i, j+1:end]) + sum(mat[j+1:end,1:i])
		end
	end
	n = diag(N)   # Extract the diagonal elements as a vector
	D = diagm(n)  # Construct a diagnonal matrix from a vector
	R = zeros(Int64, r-1)
	for i = 1:r-1
		R[i] = sum(mat[1:i,i+1:end])
	end
    #If N is a singular matrix, return 1.0 directly
	if abs(det(N)) <= eps()
        return 1.0, 0, 0, 0, 0
    else
        Ni = inv(N)
		ω2 = Ni*n   
		nu = 1.0/(n'*ω2)        # p.451
		ω1 = (D*ω2) .*nu
        Δ1 = ω1'*log.((R .+ 0.5)./(n - R .+ 0.5))
        Δ2 = log((0.5 .+ ω2' * R)/(0.5 .+ ω2' * (n .- R)))
        v1 = 4 * (1 + 0.25 * Δ1^2)*nu
        v2 = 4 * (1 + 0.25 * Δ2^2)*nu
        se = sqrt((v1 + v2) * 0.5)
        z1 = Δ1/se
		pval = pvalue(Normal(0, 1), z1, tail = :both)
		#println(ω1, ω2, Δ1, Δ2, se, z1, pval)
        return pval, Δ1, Δ2, se, z1
    end
end

mat = [43 8 3 0; 2 2 5 3; 1 0 7 2; 0 0 1 5]
McCullagh_test(mat)