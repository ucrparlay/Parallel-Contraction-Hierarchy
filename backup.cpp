
template <class Func>
auto iterate(NodeId u, Func &&f, bool forward = true)
{
    const &offset = forward ? G.offset : G.backward_offset;
    const &E = foward ? G.E : G.backward_E;
    const &edges_hash_map = foward ? out_edges_hash_map : in_edges_hash_map;
    for (size_t i = offset[u]; i < offset[u + 1], i++)
    {
        NodeId v = E[i].v;
        EdgeTy w = E[i].w;
        f(u, v, w);
    }

    // TODO: distinguish forward/backward
    unsigned long index = edges_hash_map.first_index(u);
    while (index)
    {
        if (edges_hash_map.H[index].valid)
        {
            if (edges_hash_map.H[index].key == UINT_N_MAX)
                return true;
            if (edges_hash_map.H[index].key == u)
            {
                NodeId v = edges_hash_map.H[index].value.first;
                EdgeTy w = edges_hash_map.H[index].value.second;
                f(u, v, w);
            }
        }
        index = shortcutHash.next_index(index);
    }
}