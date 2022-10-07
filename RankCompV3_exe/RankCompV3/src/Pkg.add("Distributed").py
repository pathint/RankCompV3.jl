Pkg.add("Distributed")
Pkg.add("SharedArrays")
Pkg.add("MultipleTesting")
Pkg.add("DataFrames")
Pkg.add("DelimitedFiles")
Pkg.add("CSV")
Pkg.add("Statistics")
Pkg.add("RCall")
Pkg.add("HypothesisTests")
Pkg.add("Distributions")
Pkg.add("Random")
Pkg.add("LinearAlgebra")
Pkg.add("ArgParse")

print("hello")



module Example
export hello, domath

"""
    hello(who::String)
Return "Hello, `who`".
"""
hello(who::String) = "Hello, $who"

"""
    domath(x::Number)
Return `x + 5`.
"""
domath(x::Number) = x + 5

end