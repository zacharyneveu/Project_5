
using LightGraphs
using MetaGraphs
using GraphPlot
using Glob

function loadGraph(fn)
    open(fn) do f
        s = read(f, String)
        s = split(s, '\n')
        nCols = parse(Int, s[1])
        nv = parse(Int, split(s[2], ' ')[1])
        g = SimpleGraph(nv)
        for i=3:length(s)-1
            start = parse(Int, split(s[i], ' ')[1])+1
            fin =   parse(Int, split(s[i], ' ')[2])+1
            add_edge!(g, start, fin)
        end
        mg = MetaGraph(g)
        set_prop!(mg, :nColors, nCols)
        set_prop!(mg, :name, split(fn, '/')[end])
        return mg
    end
end

function writeSolution(g, folder)
    n = get_prop(g, :name)
    n = n[1:end-6]*".output"
    nc = get_prop(g, :nColors)
    open(folder*n, "w+") do f
        confs = getConflicts(g)
        write(f, "Number of Colors: $(nc)\n")
        write(f, "Best solution: $(confs) conflicts\n")
        for v in vertices(g)
            write(f, "Node $(v): $(get_prop(g, v, :color))\n")
        end
    end
end

function writeImage(g, folder)
    n = get_prop(g, :name)
    n = n[1:end-6]*".output"
    nc = get_prop(g, :nColors)
    fillc = distinguishable_colors(nc, colorant"blue")
    elabs = [get_prop(g, e.src, :color) == get_prop(g, e.dst, :color) ? "x" : "" for e in edges(g)]
    cs = [get_prop(g, v, :color) for v in vertices(g)]
    c = [fillc[i] for i in cs]
    draw(PNG(folder*n, 16cm, 16cm), gplot(g, nodefillc=c, nodelabel=1:nv(g), edgelabel=elabs))
end

function init!(g)
    for v in vertices(g)
        nc = get_prop(g, :nColors)
        color = rand(0:nc-1)
        set_prop!(g, v, :color, color)
    end
end

function init(g)
    cg = copy(g)
    init!(cg)
    return cg
end

function plotColors(g)
    colors = [get_prop(g, v, :color) for v in vertices(g)]
    elabs = [get_prop(g, e.src, :color) == get_prop(g, e.dst, :color) ? "x" : "" for e in edges(g)]
    gplot(g, nodelabel=colors, edgelabel=elabs)
end

function getConflicts(g)
    return sum([get_prop(g, e.src, :color) == get_prop(g, e.dst, :color) for e in edges(g)])
end

function getBestNeighbor(g)
    best = g
    for v in vertices(g)
        new = copy(g)
        newc = (get_prop(new, v, :color) + 1) % get_prop(g, :nColors)
        set_prop!(new, v, :color, newc)
        if getConflicts(new) <= getConflicts(best)
            best = new
        end
    end
    if best == g
        return false
    else
        return best
    end
end
        

function steepestDescent(g)
    cg = copy(g)
    while true
        new = getBestNeighbor(cg)
        if new != false
            cg = new
        else
            return cg
        end
    end
end

using Glob

fs = glob("color*.input", "instances")
getNum(str) = parse(Int, str[16:end-8])
sort!(fs, lt=(a,b)->getNum(a) < getNum(b))

for f in fs
    fn = split(f, '/')[end]
    g = loadGraph(f)
    init!(g)
    g = steepestDescent(g)
    confs = getConflicts(g)
    writeSolution(g, "sdoutputs/")
    println("Instance $fn: $(confs) conflicts")
end 

function initPop(g, size)
    pop = [init(g) for i=1:size]
    return pop
end

function fitness(g)
    return ne(g) - getConflicts(g)
end

function mutate(g)
    v = rand(1:nv(g))
    nc = get_prop(g, :nColors)
    c = rand(0:nc-1)
    #println("Changing node $v to color $c")
end

function mate(g1, g2)
    cross = rand(1:nv(g1))
    child = copy(g1)
    for v=1:nv(g1)
        g1c = get_prop(g1, v, :color)
        g2c = get_prop(g2, v, :color)
        if v <= cross
            set_prop!(child, v, :color, g1c)
        else
            set_prop!(child, v, :color, g2c)
        end
    end
    return child
end

function select(pop, pct)
    fits = [fitness(g) for g in pop]
    max_fit = maximum(fits)
    min_fit = minimum(fits)
    fitr = max_fit-min_fit
    sel = [g for g in pop if rand(0:max_fit) < (fitness(g))*1. *pct/100]
    #println("Num selected: $(length(sel))")
    return sel
end

using StatsBase

function generate(pop,survive_rate,mutate_rate)
    # TODO: make this a do->while once I know if this is possible in Julia
    parents = select(pop, survive_rate)
    # make sure more than two parents present or mating fails...
    while length(parents) < 2
        parents = select(pop, survive_rate)
    end
    children = copy(pop)
    fits = [fitness(g) for g in parents]
    # new pop should be same size as old pop
    # TODO: experiment with this assumption
    for i=1:length(pop)
        ps = sample(parents, Weights(fits), 2, replace=false)
        children[i] = mate(ps[1], ps[2])
        if rand(0:100) > mutate_rate
            mutate(children[i])
        end
    end
    return children
end

using ProgressMeter

function evolve(pop, generations, survive_rate, mutate_rate)
    pc = copy(pop)
    @showprogress 0.5 "Instance progress: " for i=1:generations
        pc = generate(pc, survive_rate, mutate_rate)
        maxfit = maximum([fitness(g) for g in pc])
        if maxfit == ne(pop[1]) # perfect solution found
            println("Perfect solution found at generation $i")
            return pc
        end
    end
    return pc
end

meanFit(pop) = mean([fitness(g) for g in pop])
maxFit(pop) = maximum([fitness(g) for g in pop])

fs = glob("color*.input", "instances")
getNum(str) = parse(Int, str[16:end-8])
sort!(fs, lt=(a,b)->getNum(a) < getNum(b))

function run_all(fs)
    for f in fs
        fn = split(f, '/')[end]
        g = loadGraph(f)
        pop = initPop(g, 100)
        bestpop = evolve(pop, 400, 80, 40)
        inst = argmax([fitness(g) for g in bestpop])
        confs = getConflicts(bestpop[inst])
        writeSolution(bestpop[inst], "gaoutputs/")
        println("Instance $fn: $(confs) conflicts")
        print("")
    end 
end

run_all(fs)


