#!/usr/bin/env python3
"""Plots for the hands-on from Microbializer output.

Usage:
  python phylo_plots.py --og orthogroups_annotated.csv --ani ani_pairwise_values.csv \
         --tree final_species_tree.newick --out plots/

Makes: gcs_vs_ani.png, saturation.png, ushape.png, tree_matrix.png
Needs: numpy pandas matplotlib (pip install numpy pandas matplotlib)
"""
import argparse, os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument('--og', required=True, help='orthogroups_annotated.csv (one column per genome with gene IDs)')
ap.add_argument('--ani', help='ani_pairwise_values.csv from Microbializer')
ap.add_argument('--tree', help='Newick species tree')
ap.add_argument('--out', default='plots')
ap.add_argument('--perms', type=int, default=100, help='random genome orders for the saturation curve')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

# ---------- presence/absence matrix (genomes x orthogroups) ----------
df = pd.read_csv(a.og, low_memory=False)
meta = {'OG_name','OG size (number of genomes)','OG size (total number of genes)','CAI_mean','knum','knum_description'}
gcols = [c for c in df.columns if c not in meta]
pa = df[gcols].notna() & (df[gcols].astype(str).apply(lambda s: s.str.strip()) != '')
pa = pa.astype(int); pa.columns = [re.match(r'(GC[AF]_\d+\.\d)', c).group(1) if re.match(r'GC[AF]_', c) else c for c in gcols]
M = pa.values.T            # genomes x families
G, F = M.shape
print(f'{G} genomes, {F} orthogroups')

# ---------- 1. gene-content similarity vs ANI ----------
inter = M @ M.T; n = M.sum(1)
union = n[:, None] + n[None, :] - inter
jac = inter / union                                   # Jaccard gene-content similarity
ovl = inter / np.minimum(n[:, None], n[None, :])      # overlap: robust to incomplete genomes
if a.ani:
    ani = pd.read_csv(a.ani, index_col=0)
    ani = ani[[c for c in ani.columns if c in pa.columns]]
    names = [g for g in pa.columns if g in ani.index and g in ani.columns]
    idx = [list(pa.columns).index(g) for g in names]
    A = (ani.loc[names, names].values + ani.loc[names, names].values.T) / 2
    iu = np.triu_indices(len(names), 1)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), sharex=True)
    for axx, S, lab in ((ax[0], jac, 'Jaccard'), (ax[1], ovl, 'Overlap')):
        s = S[np.ix_(idx, idx)][iu]
        axx.scatter(A[iu], s, s=6, alpha=.4, color='#17695D')
        axx.set_xlabel('ANI (%)'); axx.set_ylabel(f'Gene-content similarity ({lab})')
        r = np.corrcoef(A[iu], s)[0, 1]; axx.set_title(f'{lab}: Pearson r = {r:.2f}')
    plt.tight_layout(); plt.savefig(f'{a.out}/gcs_vs_ani.png', dpi=150); plt.close()
    print('gcs_vs_ani.png')

# ---------- 2. saturation (rarefaction) curves ----------
rng = np.random.default_rng(1)
pan = np.zeros((a.perms, G)); core = np.zeros((a.perms, G))
for p in range(a.perms):
    o = rng.permutation(G); cs = np.cumsum(M[o], axis=0)
    pan[p] = (cs > 0).sum(1); core[p] = (cs == np.arange(1, G + 1)[:, None]).sum(1)
x = np.arange(1, G + 1)
fig, ax = plt.subplots(figsize=(6.5, 4.6))
for arr, c, l in ((pan, '#17695D', 'Pangenome'), (core, '#A8500F', 'Core')):
    ax.plot(x, arr.mean(0), color=c, lw=2, label=l); ax.fill_between(x, arr.min(0), arr.max(0), color=c, alpha=.15)
ax.set_xlabel('Genomes added'); ax.set_ylabel('Gene families'); ax.legend()
plt.tight_layout(); plt.savefig(f'{a.out}/saturation.png', dpi=150); plt.close()
# Heaps' law exponent gamma (n ~ k N^gamma): gamma<1 => open pangenome (alpha=1-gamma<1)
g_ = np.polyfit(np.log(x[1:]), np.log(pan.mean(0)[1:]), 1)[0]
print(f"saturation.png; Heaps' gamma = {g_:.2f} (gamma > 0 and alpha = 1-gamma < 1: open pangenome)")

# ---------- 3. U-shape ----------
cnt = M.sum(0); h = np.bincount(cnt, minlength=G + 1)[1:]
fig, ax = plt.subplots(figsize=(7, 4.4))
ax.bar(np.arange(1, G + 1), h, color='#17695D'); ax.set_yscale('log')
ax.set_xlabel('Number of genomes carrying the family'); ax.set_ylabel('Gene families (log)')
plt.tight_layout(); plt.savefig(f'{a.out}/ushape.png', dpi=150); plt.close(); print('ushape.png')

# ---------- 4. tree next to the matrix ----------
if a.tree:
    t = open(a.tree).read().strip().rstrip(';'); pos = 0
    def parse():
        global pos
        if t[pos] == '(':
            pos += 1; ch = [parse()]
            while t[pos] == ',': pos += 1; ch.append(parse())
            pos += 1
            m = re.match(r'[^:,()]*', t[pos:]); pos += m.end()      # skip support label
            bl = 0.0
            if pos < len(t) and t[pos] == ':':
                m = re.match(r':([\d.eE-]+)', t[pos:]); bl = float(m.group(1)); pos += m.end()
            return {'ch': ch, 'bl': bl}
        m = re.match(r'([^:,()]+)(?::([\d.eE-]+))?', t[pos:]); pos += m.end()
        return {'name': m.group(1), 'bl': float(m.group(2) or 0), 'ch': []}
    root = parse(); leaves = []
    def walk(n):
        if not n['ch']: leaves.append(n['name'])
        for c in n['ch']: walk(c)
    walk(root); yy = {l: i for i, l in enumerate(leaves)}
    def coords(n, x0):
        n['x'] = x0 + n['bl']
        if not n['ch']: n['y'] = yy[n['name']]
        else:
            for c in n['ch']: coords(c, n['x'])
            n['y'] = np.mean([c['y'] for c in n['ch']])
    coords(root, 0)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 7), gridspec_kw={'width_ratios': [1, 3.2], 'wspace': .01})
    def draw(n):
        if n['ch']:
            ys = [c['y'] for c in n['ch']]; a1.plot([n['x']] * 2, [min(ys), max(ys)], 'k', lw=1)
            for c in n['ch']: a1.plot([n['x'], c['x']], [c['y']] * 2, 'k', lw=1); draw(c)
    draw(root); a1.set_ylim(len(leaves) - .5, -.5); a1.axis('off')
    order = [l for l in leaves if l in pa.columns]
    sub = pa[order]; keep = (sub.sum(axis=1) >= 2) & (sub.sum(axis=1) <= len(order) - 1)
    S = sub[keep]; S = S.iloc[np.argsort(-S.sum(axis=1).values, kind='stable')]
    a2.imshow(S.values.T, aspect='auto', cmap=matplotlib.colors.ListedColormap(['#ECEFEA', '#17695D']), interpolation='nearest')
    a2.set_yticks([]); a2.set_xticks([]); a2.set_xlabel(f'{len(S)} variable gene families (sorted by frequency)')
    plt.savefig(f'{a.out}/tree_matrix.png', dpi=130, bbox_inches='tight'); plt.close(); print('tree_matrix.png')
