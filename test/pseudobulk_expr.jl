using DataFrames

function pseudobulk_expr(group_expr::DataFrame,
						ncell_pseudo::Int,
						g_name::String7)
	r_group, c_group = size(group_expr)
	it_group = collect(Iterators.partition(sample(1:c_group , c_group , replace = false), ncell_pseudo)) # Random-shuffle, then partition
	group_expr = reduce(hcat, [sum.(eachrow( group_expr[:, i])) for i in it_group ]) # Matrix r_group x ncell_pseudo
	group_expr = DataFrame(group_expr,:auto)
	rename!(group_expr,string.(g_name,"_",names(group_expr)))
	return group_expr
end

pseudobulk_expr(DataFrame(rand(0:32, 10, 6),:auto), 3, String7("group1"))