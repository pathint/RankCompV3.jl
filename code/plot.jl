using DataFrames
using StatsPlots
using Plots.PlotMeasures

function plot_result(result::DataFrame, 
		fn_stem::AbstractString = "RankCompV3")
	gr(size=(400,300))
	p1 = @df result density([:padj, :pval],
							trim = true,
							xlabel = "P-value",
							ylabel = "Density",
							label  = ["FDR" "P-Value"],
							title = "Distribution of P-values",
							titlefontsize = 10,
							fontfamily="Arial",
							xlim  = (0, 1),
							grid = false,
							margin = 8mm,
							framestype = :box,
							linewidth = 2
							);
	savefig(p1, join([fn_stem, "p_value.pdf"], "_"))

	gr(size=(900,600))
	p2_11 = histogram(result.n11, bin = 100, title = "n11");
	p2_12 = histogram(result.n12, bin = 100, title = "n12");
	p2_13 = histogram(result.n13, bin = 100, title = "n13");
	p2_21 = histogram(result.n21, bin = 100, title = "n21");
	p2_22 = histogram(result.n22, bin = 100, title = "n22");
	p2_23 = histogram(result.n23, bin = 100, title = "n23");
	p2_31 = histogram(result.n31, bin = 100, title = "n31");
	p2_32 = histogram(result.n32, bin = 100, title = "n32");
	p2_33 = histogram(result.n33, bin = 100, title = "n33");
	p2 = plot(p2_11, p2_12, p2_13, p2_21, p2_22, p2_23, p2_31,p2_32, p2_33,
			  layout = (3, 3),
			  grid = false,
			  titlefontsize = 10,
			  tickfontsize  = 7,
			  fontfamily="Arial",
			  legend = false,
			  left_margin = 12mm,
			  framestyle = :box,
			  linewidth = 0);
	savefig(p2, join([fn_stem, "contigency_table.pdf"], "_"))

	gr(size=(400,300))
	p3 = @df result density([:Δ1, :Δ2],
							trim = true,
							xlabel = "delta",
							ylabel = "Density",
							label  = ["delta1" "delta2"],
							title = "Distribution of delta-values",
							titlefontsize = 10,
							#fontfamily="Arial",
							grid = false,
							margin = 8mm,
							framestyle = :box,
							linewidth = 2
							);
	savefig(p3, join([fn_stem, "delta_value.pdf"], "_"))

	p4 = @df result density([:se],
							trim = true,
							xlabel = "s.e.",
							ylabel = "Density",
							#label  = ["s.e."],
							title = "Distribution of standard errors",
							titlefontsize = 10,
							#fontfamily="Arial",
							grid = false,
							margin = 8mm,
							framestyle = :box,
							linewidth = 2
							);
	savefig(p4, join([fn_stem, "se.pdf"], "_"))

	p5 = @df result density([:z1],
							trim = true,
							xlabel = "z",
							ylabel = "Density",
							#label  = ["z1"],
							title = "Distribution of z1",
							titlefontsize = 10,
							#fontfamily="Arial",
							grid = false,
							margin = 8mm,
							framestyle = :box,
							linewidth = 2
							);
	savefig(p5, join([fn_stem, "z1.pdf"], "_"))
	return p1, p2, p3, p4, p5
end

function plot_heatmap(df_ctrl::DataFrame, 
		df_treat::DataFrame, 
		fn_stem::AbstractString = "RankCompV3";
		  log1p::Bool = true)
	if log1p
		ctrl  = log10.(1 .+ Matrix( df_ctrl[:, Not(1)] ))
		treat = log10.(1 .+ Matrix(df_treat[:, Not(1)]))
	else
		ctrl  = Matrix( df_ctrl[:, Not(1)])
		treat = Matrix(df_treat[:, Not(1)])
	end
	gr(size=(900,600))
	phc = heatmap(ctrl,  title = "Ctrl",  yticks = false)
	pht = heatmap(treat, title = "Treat", yticks = false)
	ph = plot(phc, pht,
			  layout = (1, 2),
			  grid = false,
			  titlefontsize = 10,
			  tickfontsize  = 7,
			  #fontfamily="Arial",
			  #legend = false,
			  left_margin = 12mm,
			  linewidth = 0);

	savefig(ph, join([fn_stem, "expr_heat.pdf"], "_"))

	gr(size=(400,300))
	pd = density(hcat(ctrl, treat),
							trim = true,
							xlabel = "Expr",
							ylabel = "Density",
							legend = false,
							#label  = ["FDR" "P-Value"],
							title = "Distribution of Expression Values",
							titlefontsize = 10,
							#fontfamily="Arial",
							#xlim  = (0, 1),
							grid = false,
							margin = 8mm,
							linewidth = 2
							);
	savefig(pd, join([fn_stem, "expr_dist.pdf"], "_"))
	return ph, pd
end
