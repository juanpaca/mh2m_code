import os, sys
from collections import Counter

# =====================================================================
# Triangular Macro Mesh Generator for the Full Unit Square [0,1]^2.
#
# Each cell (i,j) of an N×N grid is split into 2 triangles.
# Each triangle is one macro element with 3 boundary edges.
#
# Level cM=0: grid 1x1 → 2 elements, h=1.0
# Level cM=1: grid 2x2 → 8 elements, h=0.5
# Level cM=2: grid 4x4 → 32 elements, h=0.25
#
# Labels: 100 + 3*K + e  (e = 0,1,2)
# =====================================================================

def generate(cM, output_dir):
    N = 2**cM                # cells per side
    h = 1.0 / N
    nn = N + 1
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

    # Each square cell → 2 triangles → 2 macro elements (nTriPerEl = 1)
    nEdgesPerEl = 3
    total_elements = 2 * N * N

    elem_data = []  # [(K, vlist_global, tris_local, bnd_local)]

    for j in range(N):
        for i in range(N):
            v00, v10, v11, v01 = (i,j), (i+1,j), (i+1,j+1), (i,j+1)

            # ---- Element K0 = diagonal (i,j)-(i+1,j)-(i+1,j+1) ----
            K0 = 2 * (j * N + i)
            verts0 = [v00, v10, v11]
            vlist0 = sorted(verts0)
            local0 = {v: idx+1 for idx, v in enumerate(vlist0)}
            tris_local0 = [tuple(local0[v] for v in verts0)]
            bnd_local0 = []
            for e in range(3):
                va, vb = verts0[e], verts0[(e+1)%3]
                label = 100 + nEdgesPerEl*K0 + e
                bnd_local0.append((local0[va], local0[vb], label))
            elem_data.append((K0, vlist0, tris_local0, bnd_local0))

            # ---- Element K1 = diagonal (i,j)-(i+1,j+1)-(i,j+1) ----
            K1 = 2 * (j * N + i) + 1
            verts1 = [v00, v11, v01]
            vlist1 = sorted(verts1)
            local1 = {v: idx+1 for idx, v in enumerate(vlist1)}
            tris_local1 = [tuple(local1[v] for v in verts1)]
            bnd_local1 = []
            for e in range(3):
                va, vb = verts1[e], verts1[(e+1)%3]
                label = 100 + nEdgesPerEl*K1 + e
                bnd_local1.append((local1[va], local1[vb], label))
            elem_data.append((K1, vlist1, tris_local1, bnd_local1))

    # ---- Global triangulation for calP.msh ---------------------------
    all_tris_global = []
    for j in range(N):
        for i in range(N):
            v00, v10, v11, v01 = gid(i,j), gid(i+1,j), gid(i+1,j+1), gid(i,j+1)
            all_tris_global.append((v00, v10, v11))
            all_tris_global.append((v00, v11, v01))

    # Boundary edges: edges appearing only once
    edge_cnt = Counter()
    for (a, b, c) in all_tris_global:
        for e in [(a,b), (b,c), (c,a)]:
            key = e if e[0] < e[1] else (e[1], e[0])
            edge_cnt[key] += 1
    bnd_global = [e for e, cnt in edge_cnt.items() if cnt == 1]

    # ---- Write calP.msh ---------------------------------------------
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

    print(f"calP: {nnodes} nodes, {len(all_tris_global)} tris, "
          f"{total_elements} tri elements, {len(bnd_global)} bnd edges")

    # ---- Write Th_K.msh files ---------------------------------------
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

    print(f"Generated {len(elem_data)} Th_K.msh files in '{output_dir}/'")

if __name__ == "__main__":
    cM = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    out = f"meshes_tri_cM{cM}"
    generate(cM, out)
