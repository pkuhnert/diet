tree.depth <-
function (nodes)
{
    depth <- floor(log(nodes, base = 2) + 1e-07)
    as.vector(depth - min(depth))
}

