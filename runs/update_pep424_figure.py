"""Rebuild the published-prediction comparison from committed result values.

The FoldAmyloid and transfer bars use the same 419-record pep424 alignment.
The AmyloGram 0.865 bar is its separate published cross-validation protocol,
not a same-split head-to-head score.
"""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[1]
result = json.loads((root / 'results/pep424_v3.json').read_text())
values = [result['foldamyloid_auc'], result['transfer_pepx_gnn_auc'],
          result['cv5_pepx_gnn_auc'], result['amylogram_published_cv_auc']]
labels = ['FoldAmyloid\n(published preds)', 'pepx-GNN transfer\n(same 419)',
          'pepx-GNN\n(5-fold CV)', 'AmyloGram\n(published CV)']
fig, ax = plt.subplots(figsize=(7.4, 4.8), dpi=160)
bars = ax.bar(labels, values, color=['#999999', '#287c5a', '#3a986b', '#999999'])
ax.set_ylim(0.5, 0.93)
ax.set_ylabel('AUROC')
ax.set_title('pep424 amyloid benchmark: aligned transfer and separate CV')
for bar, val in zip(bars, values):
    ax.annotate(f'{val:.3f}', (bar.get_x()+bar.get_width()/2, val),
                ha='center', va='bottom', xytext=(0, 4), textcoords='offset points')
ax.spines[['top','right']].set_visible(False)
fig.tight_layout()
fig.savefig(root / 'paper/figures/fig2_pep424_headtohead.png')
plt.close(fig)
