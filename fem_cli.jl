# fem_cli.jl — minimal & robust
# Usage:
#   julia fem_cli.jl IN_PARQUET FORMULA_TXT CLUSTERS METHOD OUT_PREFIX WEIGHT_VAR WEIGHT_TYPE [WARM_N]
# Example:
#   julia fem_cli.jl temp.parquet formula.txt ntrab cuda fem_out weight_var aweight 200000
#   julia fem_cli.jl temp.parquet formula.txt ""    cpu  fem_out "" ""

@time using DataFrames, StatsModels, FixedEffectModels, Vcov, LinearAlgebra, CSV, Parquet2, CategoricalArrays, CUDA, Dates
#BLAS.set_num_threads(1)

# I want to know the time of the next lines

# @time begin
#     # --- load CUDA at top-level if available ---
#     const HAVE_CUDA = try
#         @eval using CUDA
#         true
#     catch
#         false
#     end
# end

# println("Current time a: ", Dates.format(now(), "HH:MM:SS"))

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

# println("Current time b: ", Dates.format(now(), "HH:MM:SS"))

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

# println("Current time c: ", Dates.format(now(), "HH:MM:SS"))

function main(args)
    # println("Current time 1: ", Dates.format(now(), "HH:MM:SS"))

    length(args) ≥ 7 || error("Need: IN_PARQUET FORMULA_TXT CLUSTERS METHOD OUT_PREFIX WEIGHT_VAR WEIGHT_TYPE [WARM_N]")
    parquet_path, formula_txt, clusters_str, method_str, outprefix, weight_var_raw, weight_type_raw = args[1:7]
    warm_n = (length(args) ≥ 8 && tryparse(Int, args[8]) !== nothing) ? parse(Int, args[8]) : 0
    weight_var  = strip(weight_var_raw)
    weight_type = lowercase(strip(weight_type_raw))
    if isempty(weight_var) != isempty(weight_type)
        @warn "Weight variable/type mismatch; ignoring weights." weight_var weight_type
        weight_var = ""
        weight_type = ""
    end

    # println("Current time 2: ", Dates.format(now(), "HH:MM:SS"))

    # read + translate formula
    raw_formula = strip(read(formula_txt, String))
    translated, ivars, ibase, fevars = translate_stata_factors(raw_formula)
    println("Formula (Stata):        ", raw_formula)
    println("Formula (StatsModels):  ", translated)
    println("Clusters:               ", clusters_str)
    if isempty(weight_type) || isempty(weight_var)
        println("Weights:                (none)")
    else
        println("Weights:                ", weight_type, " -> ", weight_var)
    end

    # println("Current time 3: ", Dates.format(now(), "HH:MM:SS"))

    # load only referenced columns (RHS+LHS+FE+clusters)
    cols = Symbol.(needed_vars(translated, fevars, clusters_str))
    weight_sym = nothing
    if !isempty(weight_var)
        try
            weight_sym = Symbol(weight_var)
            if !(weight_sym in cols)
                push!(cols, weight_sym)
            end
        catch err
            error("Invalid weight variable name '$weight_var'")
        end
    end
    @time ds = Parquet2.Dataset(parquet_path)
    @time df = DataFrame(Parquet2.select(ds, cols...); copycols=false)
    if weight_sym !== nothing
        if !hasproperty(df, weight_sym)
            error("Weight column '$weight_var' not found in parquet data")
        end
        wcol_any = df[!, weight_sym]
        if !(eltype(wcol_any) <: Real)
            error("Weight column '$weight_var' must be numeric")
        end
        wvals = Float64.(wcol_any)
        if weight_type in ("fweight", "frequency")
            if any(x -> !isfinite(x) || x < 0 || x != round(x), wvals)
                error("fweights must be nonnegative integers")
            end
            df[!, weight_sym] = round.(Int, wvals)
        elseif weight_type in ("pweight", "probability")
            if any(x -> !isfinite(x) || x <= 0, wvals)
                error("pweights must be positive")
            end
            df[!, weight_sym] = wvals
        else
            if any(x -> !isfinite(x) || x <= 0, wvals)
                @warn "Non-positive analytic weights detected; they will be passed through as-is." weight_var
            end
            df[!, weight_sym] = wvals
        end
    end
    # display the column names loaded
    println("Loaded columns: ", names(df))

    # println("Current time 4: ", Dates.format(now(), "HH:MM:SS"))

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

    # println("Current time 5: ", Dates.format(now(), "HH:MM:SS"))

    # build formula
    f = eval(Meta.parse("@formula($(translated))"))

    # VCOV: robust if clusters empty/keyword, else cluster(...)
    clusters = strip(clusters_str)
    vc = if isempty(clusters) || lowercase(clusters) in ("robust","none","-","0","false")
        Vcov.robust()
    else
        Vcov.cluster(Symbol.(strip.(split(clusters, ",")))...)
    end

    # println("Current time 6: ", Dates.format(now(), "HH:MM:SS"))

    # method mapping
    ms = lowercase(strip(method_str))
    method = ms in ("cuda","gpu") ? :CUDA :
             ms == "cpu"          ? :cpu  :
             error("METHOD must be CPU or CUDA/GPU, got '$method_str'")

    # optional GPU warm-up
    if method === :CUDA
        #HAVE_CUDA || error("METHOD=CUDA requested but CUDA.jl is not available")
        @assert CUDA.functional() "CUDA is not functional on this machine"
        CUDA.precompile_runtime()
        if warm_n > 0
            warm = first(df, min(nrow(df), warm_n))
            warm_kwargs = weight_sym === nothing ? (; progress_bar=false) : (; weights=weight_sym, progress_bar=false)
            @time _ = reg(warm, f; method=:CUDA, warm_kwargs...)
            CUDA.synchronize()
        end
        println("CUDA is ready.")
        println("CPU threads: ", Threads.nthreads())
        println("BLAS threads: ", BLAS.get_num_threads())
    end

    # println("Current time 7: ", Dates.format(now(), "HH:MM:SS"))

    # fit + outputs
    reg_kwargs = weight_sym === nothing ? (; progress_bar=false) : (; weights=weight_sym, progress_bar=false)
    @time m = reg(df, f, vc; method=method, reg_kwargs...)
    open(outprefix * "_summary.txt", "w") do io
        show(io, MIME"text/plain"(), m); println(io)
    end

    # println("Current time 8: ", Dates.format(now(), "HH:MM:SS"))

    β = coef(m); V = Matrix(vcov(m)); se = sqrt.(diag(V))
    CSV.write(outprefix * "_coef.csv", DataFrame(term=StatsModels.coefnames(m), estimate=β, stderr=se))

    stats_rows = NamedTuple{(:stat, :value), Tuple{String, Float64}}[]
    function push_stat!(label::AbstractString, thunk::Function)
        val = try
            thunk()
        catch err
            @warn "Could not compute $label" exception=(err, catch_backtrace())
            nothing
        end
        if val isa Number && isfinite(val)
            push!(stats_rows, (stat=String(label), value=float(val)))
        end
        return nothing
    end

    push_stat!("r2", () -> FixedEffectModels.r2(m))
    if isdefined(FixedEffectModels, :adjr2)
        push_stat!("adj_r2", () -> FixedEffectModels.adjr2(m))
    elseif isdefined(StatsModels, :adjr2)
        push_stat!("adj_r2", () -> StatsModels.adjr2(m))
    end
    if isdefined(FixedEffectModels, :r2_within)
        push_stat!("r2_within", () -> FixedEffectModels.r2_within(m))
    end
    if isdefined(FixedEffectModels, :r2_between)
        push_stat!("r2_between", () -> FixedEffectModels.r2_between(m))
    end

    if !isempty(stats_rows)
        CSV.write(outprefix * "_stats.csv", DataFrame(stats_rows))
    end

    # println("Current time 9: ", Dates.format(now(), "HH:MM:SS"))
end

main(ARGS)
