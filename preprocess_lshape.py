import os, sys
from collections import Counter

# =====================================================================
# L-Shaped Macro Mesh Generator for the Full Unit Square [0,1]^2.
#
# Base level (cM=0):  4 elements (2^2),  grid 4x4 cells,  h=0.25
# Level cM=1:         16 elements (2^4), grid 8x8 cells,  h=0.125
# Level cM=2:         64 elements (2^6), grid 16x16 cells, h=0.0625
#
# Element layout per 4x4 macro-block (cell offsets from block corner):
#
#   gi=0 (bottom-left L):  cells (0,0),(1,0),(1,1),(1,2)
#   gi=1 (top-left L):     cells (0,1),(0,2),(0,3),(1,3)
#   gi=2 (top-right L):    cells (2,3),(3,1),(3,2),(3,3)
#   gi=3 (bottom-right L): cells (2,0),(2,1),(2,2),(3,0)
#
# Each element K has 6 boundary edges (hexagon) with labels:
#   100 + 6*K + e   (e = 0..5)
# =====================================================================

# ---- Element type definitions ----------------------------------------
# For each gi: 4 cell offsets (ci,cj) and 6 boundary edge vertex-pair
# offsets (di,dj) in CCW order.

TYPE_CELLS = [
    [(0,0), (1,0), (1,1), (1,2)],
    [(0,1), (0,2), (0,3), (1,3)],
    [(2,3), (3,1), (3,2), (3,3)],
    [(2,0), (2,1), (2,2), (3,0)],
]

TYPE_BND = [
    # gi=0
    [((0,0),(2,0)), ((2,0),(2,3)), ((2,3),(1,3)),
     ((1,3),(1,1)), ((1,1),(0,1)), ((0,1),(0,0))],
    # gi=1
    [((0,1),(1,1)), ((1,1),(1,3)), ((1,3),(2,3)),
     ((2,3),(2,4)), ((2,4),(0,4)), ((0,4),(0,1))],
    # gi=2
    [((4,4),(2,4)), ((2,4),(2,3)), ((2,3),(3,3)),
     ((3,3),(3,1)), ((3,1),(4,1)), ((4,1),(4,4))],
    # gi=3
    [((4,1),(3,1)), ((3,1),(3,3)), ((3,3),(2,3)),
     ((2,3),(2,0)), ((2,0),(4,0)), ((4,0),(4,1))],
]

def generate(cM, output_dir):
    nb = 4 * (2**cM)       # cells per side
    h = 1.0 / nb
    nn = nb + 1
    os.makedirs(output_dir, exist_ok=True)

    # ---- Build global vertex grid ------------------------------------
    node_xyz = {}
    node_id = {}
    nid = 0
    for j in range(nn):
        for i in range(nn):
            nid += 1
            node_xyz[nid] = (i * h, j * h, 0.0)
            node_id[(i, j)] = nid
    nnodes = nid

    def gid(i, j):
        return node_id[(i, j)]

    nbx = nb // 4           # macro-blocks per side
    nTriPerEl = 8           # 4 cells x 2 triangles
    total_elements = 4 * nbx * nbx

    # Collect per-element data for calP and Th_K
    all_tris_global = []    # [[(a,b,c),...] for each element] — gid-based
    elem_data = []          # [(K, verts_global, tris_local, bnd_local)] per element

    for by in range(nbx):
        for bx in range(nbx):
            i0 = 4 * bx
            j0 = 4 * by

            for gi in range(4):
                K = gi + 4 * (bx + nbx * by)
                cells = TYPE_CELLS[gi]
                bnd_pairs = TYPE_BND[gi]

                # ---- Collect all (global cell) vertices for this element ----
                vert_global_cell = set()  # set of (gi,gj)
                for (ci, cj) in cells:
                    vert_global_cell.add((i0 + ci,     j0 + cj))
                    vert_global_cell.add((i0 + ci + 1, j0 + cj))
                    vert_global_cell.add((i0 + ci + 1, j0 + cj + 1))
                    vert_global_cell.add((i0 + ci,     j0 + cj + 1))

                # ---- Triangles (global cell coords) ----
                tris_cell = []  # list of ((gi,gj), (gi,gj), (gi,gj))
                for (ci, cj) in cells:
                    v00 = (i0 + ci,     j0 + cj)
                    v10 = (i0 + ci + 1, j0 + cj)
                    v11 = (i0 + ci + 1, j0 + cj + 1)
                    v01 = (i0 + ci,     j0 + cj + 1)
                    tris_cell.append((v00, v10, v11))
                    tris_cell.append((v00, v11, v01))

                # ---- Boundary edges (global cell coords) ----
                bnd_cell = []
                for (fa, fb) in bnd_pairs:
                    va = (i0 + fa[0], j0 + fa[1])
                    vb = (i0 + fb[0], j0 + fb[1])
                    bnd_cell.append((va, vb))

                # ---- Local node numbering for Th_K.msh ----
                vlist = sorted(vert_global_cell)  # deterministic order
                local_of = {v: idx+1 for idx, v in enumerate(vlist)}

                # Local triangles
                tris_local = []
                for tri in tris_cell:
                    tris_local.append(tuple(local_of[v] for v in tri))

                # Local boundary edges
                bnd_local = []
                for eidx, (va, vb) in enumerate(bnd_cell):
                    label = 100 + 6*K + eidx
                    bnd_local.append((local_of[va], local_of[vb], label))

                # Global triangles (for calP.msh)
                tris_global = [tuple(gid(v[0], v[1]) for v in tri) for tri in tris_cell]

                all_tris_global.extend(tris_global)
                elem_data.append((K, vlist, tris_local, bnd_local))

    # ---- Write calP.msh ----------------------------------------------
    # Global boundary edges: edges that appear only once in the global triangulation
    edge_cnt = Counter()
    for (a, b, c) in all_tris_global:
        for e in [(a,b), (b,c), (c,a)]:
            key = e if e[0] < e[1] else (e[1], e[0])
            edge_cnt[key] += 1
    bnd_global = [e for e, cnt in edge_cnt.items() if cnt == 1]

    with open(f"{output_dir}/calP.msh", "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write(f"$Nodes\n{nnodes}\n")
        for n in sorted(node_xyz.keys()):
            x, y, z = node_xyz[n]
            f.write(f"{n} {x:.15e} {y:.15e} {z:.15e}\n")
        f.write("$EndNodes\n")
        ne = len(all_tris_global) + len(bnd_global)
        f.write(f"$Elements\n{ne}\n")
        eid = 0
        for (a, b, c) in all_tris_global:
            eid += 1
            f.write(f"{eid} 2 2 0 0 {a} {b} {c}\n")
        for (a, b) in bnd_global:
            eid += 1
            f.write(f"{eid} 1 2 1 0 {a} {b}\n")
        f.write("$EndElements\n")

    n_actual = len(elem_data)
    print(f"calP: {nnodes} nodes, {len(all_tris_global)} tris, {n_actual} L-elements (target {total_elements}), {len(bnd_global)} bnd edges")

    # ---- Write Th_K.msh ----------------------------------------------
    for (K, vlist, tris_local, bnd_local) in elem_data:
        with open(f"{output_dir}/Th_{K}.msh", "w") as f:
            f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
            nv = len(vlist)
            f.write(f"$Nodes\n{nv}\n")
            for lid, (gi, gj) in enumerate(vlist, 1):
                x = gi * h
                y = gj * h
                f.write(f"{lid} {x:.15e} {y:.15e} 0.0\n")
            f.write("$EndNodes\n")
            le = len(tris_local) + len(bnd_local)
            f.write(f"$Elements\n{le}\n")
            eid = 0
            for (a, b, c) in tris_local:
                eid += 1
                f.write(f"{eid} 2 2 0 0 {a} {b} {c}\n")
            for (a, b, lbl) in bnd_local:
                eid += 1
                f.write(f"{eid} 1 2 {lbl} 0 {a} {b}\n")
            f.write("$EndElements\n")

    print(f"Generated {n_actual} Th_K.msh files in '{output_dir}/'")

if __name__ == "__main__":
    cM = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    out = f"meshes_lshape_cM{cM}"
    generate(cM, out)
