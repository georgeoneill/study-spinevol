function faces = tet2face(tet)

all_faces = [tet(:, [1 2 3]);
    tet(:, [1 2 4]);
    tet(:, [1 3 4]);
    tet(:, [2 3 4])];
sorted_faces = sort(all_faces, 2);
[faces, ~, ~] = unique(sorted_faces, 'rows', 'stable');