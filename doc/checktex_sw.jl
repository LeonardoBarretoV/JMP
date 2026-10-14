function check(p)
    s = read(p, String)
    problems = 0

    lines = split(s, '\n'); envs = String[]
    for (k, ln) in enumerate(lines)
        startswith(strip(ln), "%") && continue
        for m in eachmatch(r"\\begin\{([A-Za-z*]+)\}", ln); push!(envs, m.captures[1]); end
        for m in eachmatch(r"\\end\{([A-Za-z*]+)\}", ln)
            e = m.captures[1]
            if isempty(envs)
                println("line $k: \\end{$e} with nothing open"); problems += 1
            elseif envs[end] != e
                println("line $k: \\end{$e} closes \\begin{$(envs[end])}"); problems += 1
                pop!(envs)
            else
                pop!(envs)
            end
        end
    end
    if !isempty(envs); println("unclosed environments: ", envs); problems += 1; end
    println("environments: balanced")

    nb = count("{", s); ne = count("}", s)
    println("braces:  { = $nb   } = $ne   diff = $(nb-ne)")
    nb == ne || (problems += 1)

    nd = count("\$", s) - 2*count(raw"\$", s)
    println("unescaped \$ = $nd  ", iseven(nd) ? "(even, OK)" : "(ODD -- unbalanced)")
    iseven(nd) || (problems += 1)

    # package check: match the loaded lists, not exact-brace names
    loaded = String[]
    for m in eachmatch(r"\\usepackage(?:\[[^\]]*\])?\{([^}]+)\}", s)
        append!(loaded, strip.(split(m.captures[1], ",")))
    end
    println("packages loaded: ", join(loaded, ", "))
    for (cmd, pkg) in ((raw"\mathds","dsfont"), (raw"\toprule","booktabs"),
                       (raw"\begin{shaded}","framed"), (raw"\boxed","amsmath"),
                       (raw"\definecolor","xcolor"), (raw"\includegraphics","graphicx"),
                       (raw"\begin{empheq}","empheq"))
        if occursin(cmd, s) && !(pkg in loaded)
            println("MISSING \\usepackage{$pkg} for $cmd"); problems += 1
        end
    end

    defs = Set(m.captures[1] for m in eachmatch(r"\\label\{([^}]+)\}", s))
    refs = Set(m.captures[1] for m in eachmatch(r"\\(?:eq)?ref\{([^}]+)\}", s))
    miss = setdiff(refs, defs)
    if !isempty(miss); println("UNDEFINED labels referenced: ", collect(miss)); problems += 1; end
    dup = Dict{String,Int}()
    for m in eachmatch(r"\\label\{([^}]+)\}", s); dup[m.captures[1]] = get(dup,m.captures[1],0)+1; end
    rep = [k for (k,v) in dup if v > 1]
    if !isempty(rep); println("DUPLICATE labels: ", rep); problems += 1; end
    println("labels: $(length(defs)) defined, $(length(refs)) referenced")

    used = Set(m.captures[1] for m in eachmatch(r"\\textcite\{([^}]+)\}", s))
    shim = Set(m.captures[1] for m in eachmatch(r"\\defcite\{([^}]+)\}", s))
    un = setdiff(used, shim)
    isempty(un) || println("keys with no standalone text: ", collect(un))
    println("citations: $(length(used)) keys used, $(length(shim)) shimmed")

    println(problems == 0 ? "==> STRUCTURAL CHECKS PASSED" : "==> $problems PROBLEM(S)")
    return problems
end
exit(check("JMP_inflation_model.tex") == 0 ? 0 : 1)
