index_singularity = 1;
node_indexes = mod(index_singularity - 1 + [0 1 2 3 4 5],6) + 1
found = find(ismember(node_indexes,[1,2,3]))
vertex_indexes = node_indexes(found)
