# fem_cli.jl — minimal & robust
# Usage:
#   julia fem_cli.jl IN_PARQUET FORMULA_TXT CLUSTERS METHOD OUT_PREFIX [WARM_N]
# Example:
#   julia fem_cli.jl temp.parquet formula.txt ntrab cuda fem_out 200000
#   julia fem_cli.jl temp.parquet formula.txt ""    cpu  fem_out

using DataFrames, StatsModels, FixedEffectModels, Vcov, LinearAlgebra, CSV, Parquet2, CategoricalArrays
BLAS.set_num_threads(1)

# --- load CUDA at top-level if available ---
const HAVE_CUDA = try
    @eval using CUDA
    true
catch
    false
end

function translate_stata_factors(s::AbstractString)
    t = String(s)

    # 1) Protect fe(...) and collect FE vars
    fe_blocks = Dict{String,String}()
    fevars = String[]
    i = 0
    t2 = t
    for m in eachmatch(r"fe\(([^)]+)\)", t)
        i += 1
        key = "__FE_BLOCK_$(i)__"
        fe_blocks[key] = m.match           # full "fe(...)" text
        t2 = replace(t2, m.match => key)   # hide it
        push!(fevars, strip(m.captures[1]))
    end

    # 2) Collect i.vars and ib#.vars (base levels)
    ivars = String[]
    ibase = Dict{String,String}()

    # ib(3).var or ib3.var, var must start with a letter (avoid placeholders)
    for m in eachmatch(r"\bib(?:\(([^)]+)\)|(\d+))\.([A-Za-z]\w*)", t2)
        base = something(m.captures[1], m.captures[2])
        var  = m.captures[3]
        ibase[var] = strip(base)
        push!(ivars, var)
    end
    # plain i.var (var must start with a letter)
    for m in eachmatch(r"\bi\.([A-Za-z]\w*)", t2)
        push!(ivars, m.captures[1])
    end
    ivars = unique(ivars)

    # 3) Operators: Stata -> StatsModels
    t2 = replace(t2, "##" => "*")   # main + interaction
    t2 = replace(t2, "#"  => "&")   # interaction only

    # 4) Strip prefixes on *real* variable names
    t2 = replace(t2, r"\b(i|ib(?:\([^)]+\)|\d+))\.([A-Za-z]\w*)" => s"\2")
    t2 = replace(t2, r"\bc\.([A-Za-z]\w*)" => s"\1")

    # 5) Restore fe(...) blocks
    for (k, v) in fe_blocks
        t2 = replace(t2, k => v)
    end

    # 6) Sanity cleanup: remove any stray 'i.' or 'ib... .' that might precede placeholders/FE
    #    (e.g., turns 'i.__FE_BLOCK_1__' into '__FE_BLOCK_1__', which then is fe(...))
    t2 = replace(t2, r"\bib(?:\([^)]+\)|\d+)\." => "")
    t2 = replace(t2, r"\bi\." => "")

    return t2, ivars, ibase, unique(fevars)
end

# --- pull all variables we need to load from parquet (formula + FE + clusters) ---
function needed_vars(translated::String, fevars::Vector{String}, clusters_str::String)
    # remove fe(...) blocks; we already collected fevars
    t = replace(translated, r"fe\([^)]+\)" => "")
    # replace operators with spaces and extract identifiers
    t = replace(t, r"[~\+\*\&\(\)]" => " ")
    toks = [m.match for m in eachmatch(r"\b[A-Za-z_]\w*\b", t)]
    cls  = [strip(c) for c in split(clusters_str, ",") if !isempty(strip(c))]
    return unique(vcat(toks, fevars, cls))
end

function main(args)
    length(args) ≥ 5 || error("Need: IN_PARQUET FORMULA_TXT CLUSTERS METHOD OUT_PREFIX [WARM_N]")
    parquet_path, formula_txt, clusters_str, method_str, outprefix = args[1:5]
    warm_n = (length(args) ≥ 6 && tryparse(Int, args[6]) !== nothing) ? parse(Int, args[6]) : 0

    # read + translate formula
    raw_formula = strip(read(formula_txt, String))
    translated, ivars, ibase, fevars = translate_stata_factors(raw_formula)
    println("Formula (Stata):        ", raw_formula)
    println("Formula (StatsModels):  ", translated)
    println("Clusters:               ", clusters_str)

    # load only referenced columns (RHS+LHS+FE+clusters)
    cols = Symbol.(needed_vars(translated, fevars, clusters_str))
    @time ds = Parquet2.Dataset(parquet_path)
    @time df = DataFrame(Parquet2.select(ds, cols...); copycols=false)


    # make i.-vars categorical and apply ib(#) base (by string match → robust for any level type)
    for v in ivars
        s = Symbol(v)
        if !hasproperty(df, s)
            @warn "Factor variable '$v' not found in data"; continue
        end
        df[!, s] = categorical(df[!, s]; ordered=false)
        if haskey(ibase, v)
            base_str = ibase[v]
            levs = levels(df[!, s])
            if isempty(levs)
                @warn "No levels for '$v'; cannot set base"; continue
            end
            strlevs = string.(levs)
            idx = findfirst(==(base_str), strlevs)
            if idx === nothing
                let n = min(length(levs), 10)
                    @warn "Requested base $v=$base_str not found; sample: $(first(strlevs, n))"
                end
            else
                newlevs = vcat(levs[idx], levs[setdiff(1:length(levs), idx)])
                levels!(df[!, s], newlevs)
            end
        end
    end

    # build formula
    f = eval(Meta.parse("@formula($(translated))"))

    # VCOV: robust if clusters empty/keyword, else cluster(...)
    clusters = strip(clusters_str)
    vc = if isempty(clusters) || lowercase(clusters) in ("robust","none","-","0","false")
        Vcov.robust()
    else
        Vcov.cluster(Symbol.(strip.(split(clusters, ",")))...)
    end

    # method mapping
    ms = lowercase(strip(method_str))
    method = ms in ("cuda","gpu") ? :CUDA :
             ms == "cpu"          ? :cpu  :
             error("METHOD must be CPU or CUDA/GPU, got '$method_str'")

    # optional GPU warm-up
    if method === :CUDA
        HAVE_CUDA || error("METHOD=CUDA requested but CUDA.jl is not available")
        @assert CUDA.functional() "CUDA is not functional on this machine"
        CUDA.precompile_runtime()
        if warm_n > 0
            warm = first(df, min(nrow(df), warm_n))
            @time _ = reg(warm, f; method=:CUDA)
            CUDA.synchronize()
        end
    end

    # fit + outputs
    @time m = reg(df, f, vc; method=method)
    open(outprefix * "_summary.txt", "w") do io
        show(io, MIME"text/plain"(), m); println(io)
    end
    β = coef(m); V = Matrix(vcov(m)); se = sqrt.(diag(V))
    CSV.write(outprefix * "_coef.csv", DataFrame(term=StatsModels.coefnames(m), estimate=β, stderr=se))
end

main(ARGS)