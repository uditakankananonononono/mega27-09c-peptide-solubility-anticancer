"""Deep models (PyTorch, CPU-first): PeptideCNN and PeptideGNN, plus the
multi-task TriNet that shares one encoder across the three tasks
(anticancer / solubility / aggregation)."""
from __future__ import annotations

from typing import Dict, List, Optional

import torch
import torch.nn as nn
import torch.nn.functional as F

from .alphabet import AA_STANDARD

VOCAB = len(AA_STANDARD) + 1  # + X/pad
PAD = len(AA_STANDARD)


class PeptideCNN(nn.Module):
    """Multi-kernel 1D CNN over a learned residue embedding, with an optional
    physicochemical-descriptor side head (feature fusion)."""

    def __init__(self, emb_dim: int = 32, kernels: tuple = (3, 5, 7),
                 filters: int = 64, n_desc: int = 18, dropout: float = 0.3):
        super().__init__()
        self.emb = nn.Embedding(VOCAB, emb_dim, padding_idx=PAD)
        self.convs = nn.ModuleList([
            nn.Conv1d(emb_dim, filters, k, padding=k // 2) for k in kernels])
        self.n_desc = n_desc
        fused = filters * len(kernels) + n_desc
        self.head = nn.Sequential(
            nn.Dropout(dropout), nn.Linear(fused, 64), nn.ReLU(),
            nn.Dropout(dropout), nn.Linear(64, 1))

    def forward(self, idx: torch.Tensor, desc: torch.Tensor) -> torch.Tensor:
        # idx: (B, L) long; desc: (B, n_desc) float
        x = self.emb(idx).transpose(1, 2)          # (B, E, L)
        feats = [F.relu(c(x)).amax(dim=2) for c in self.convs]  # (B, F) each
        z = torch.cat(feats + [desc], dim=1)
        return self.head(z).squeeze(-1)            # logits (B,)


class GraphLayer(nn.Module):
    """One GCN layer: H' = sigma(D^-1/2 (A+I) D^-1/2 H W)."""

    def __init__(self, in_dim: int, out_dim: int):
        super().__init__()
        self.lin = nn.Linear(in_dim, out_dim)

    def forward(self, h: torch.Tensor, adj: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        b, l, _ = h.shape
        eye = torch.eye(l, device=h.device).unsqueeze(0)
        a = adj + eye
        deg = a.sum(-1).clamp(min=1).pow(-0.5)
        d_inv = torch.diag_embed(deg)
        a_norm = d_inv @ a @ d_inv
        out = torch.bmm(a_norm, self.lin(h))
        return F.relu(out) * mask.unsqueeze(-1)


def chain_adjacency(idx: torch.Tensor, window: int = 3) -> torch.Tensor:
    """Sequence-chain adjacency with a local window (|i-j| <= window), masked
    to real (non-pad) residues."""
    b, l = idx.shape
    device = idx.device
    ar = torch.arange(l, device=device)
    near = (ar.unsqueeze(0) - ar.unsqueeze(1)).abs() <= window
    real = (idx != PAD)
    adj = near.unsqueeze(0) & real.unsqueeze(1) & real.unsqueeze(2)
    return adj.float()


class PeptideGNN(nn.Module):
    """GCN over the sequence graph (chain + local window) with descriptor
    fusion - the GNN arm required by the program spec."""

    def __init__(self, emb_dim: int = 32, hidden: int = 64, layers: int = 3,
                 n_desc: int = 18, dropout: float = 0.3, window: int = 3):
        super().__init__()
        self.emb = nn.Embedding(VOCAB, emb_dim, padding_idx=PAD)
        self.layers = nn.ModuleList([GraphLayer(emb_dim if i == 0 else hidden, hidden)
                                     for i in range(layers)])
        self.window = window
        self.head = nn.Sequential(
            nn.Dropout(dropout), nn.Linear(hidden + n_desc, 64), nn.ReLU(),
            nn.Dropout(dropout), nn.Linear(64, 1))

    def forward(self, idx: torch.Tensor, desc: torch.Tensor) -> torch.Tensor:
        mask = (idx != PAD).float()
        adj = chain_adjacency(idx, self.window)
        h = self.emb(idx)
        for layer in self.layers:
            h = layer(h, adj, mask)
        pooled = (h * mask.unsqueeze(-1)).sum(1) / mask.sum(1, keepdim=True).clamp(min=1)
        z = torch.cat([pooled, desc], dim=1)
        return self.head(z).squeeze(-1)


class TriNet(nn.Module):
    """Multi-task shared CNN encoder with three heads: acp / solubility /
    aggregation. This is the discovery machine: one embedding of 'what makes
    a peptide developable and active', scored jointly at screen time."""

    def __init__(self, emb_dim: int = 48, filters: int = 96, n_desc: int = 18,
                 dropout: float = 0.3):
        super().__init__()
        self.emb = nn.Embedding(VOCAB, emb_dim, padding_idx=PAD)
        self.convs = nn.ModuleList([nn.Conv1d(emb_dim, filters, k, padding=k // 2)
                                    for k in (3, 5, 7)])
        fused = filters * 3 + n_desc
        self.shared = nn.Sequential(nn.Dropout(dropout), nn.Linear(fused, 128), nn.ReLU())
        self.heads = nn.ModuleDict({
            t: nn.Sequential(nn.Dropout(dropout), nn.Linear(128, 1))
            for t in ("acp", "sol", "agg")})

    def encode(self, idx: torch.Tensor, desc: torch.Tensor) -> torch.Tensor:
        x = self.emb(idx).transpose(1, 2)
        feats = [F.relu(c(x)).amax(dim=2) for c in self.convs]
        return self.shared(torch.cat(feats + [desc], dim=1))

    def forward(self, idx: torch.Tensor, desc: torch.Tensor,
                task: str) -> torch.Tensor:
        return self.heads[task](self.encode(idx, desc)).squeeze(-1)

    def score_all(self, idx: torch.Tensor, desc: torch.Tensor) -> Dict[str, torch.Tensor]:
        z = self.encode(idx, desc)
        return {t: torch.sigmoid(h(z)).squeeze(-1) for t, h in self.heads.items()}


class ResConvBlock(nn.Module):
    """Dilated residual conv block: x + Conv(relu(BN(Conv(x))))."""

    def __init__(self, ch: int, dilation: int):
        super().__init__()
        self.c1 = nn.Conv1d(ch, ch, 3, padding=dilation, dilation=dilation)
        self.bn = nn.BatchNorm1d(ch)
        self.c2 = nn.Conv1d(ch, ch, 3, padding=dilation, dilation=dilation)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = self.c2(F.relu(self.bn(self.c1(x))))
        return F.relu(x + h)


class PeptideCNNv2(nn.Module):
    """Residual dilated CNN + additive attention pooling + wide feature head
    (18 physicochemical + 20 AAC + 400 DPC = 438 fused features)."""

    def __init__(self, emb_dim: int = 48, ch: int = 96, blocks: int = 3,
                 n_wide: int = 438, dropout: float = 0.35):
        super().__init__()
        self.emb = nn.Embedding(VOCAB, emb_dim, padding_idx=PAD)
        self.inp = nn.Conv1d(emb_dim, ch, 1)
        self.blocks = nn.ModuleList([ResConvBlock(ch, d) for d in (1, 2, 4)[:blocks]])
        self.attn = nn.Conv1d(ch, 1, 1)
        self.head = nn.Sequential(
            nn.Dropout(dropout), nn.Linear(2 * ch + n_wide, 128), nn.ReLU(),
            nn.Dropout(dropout), nn.Linear(128, 1))

    def forward(self, idx: torch.Tensor, wide: torch.Tensor) -> torch.Tensor:
        x = self.inp(self.emb(idx).transpose(1, 2))
        for b in self.blocks:
            x = b(x)
        a = torch.softmax(self.attn(x).masked_fill((idx == PAD).unsqueeze(1), -1e9), dim=2)
        pooled_attn = torch.bmm(a, x.transpose(1, 2)).squeeze(1)
        pooled_max = x.amax(dim=2)
        z = torch.cat([pooled_attn, pooled_max, wide], dim=1)
        return self.head(z).squeeze(-1)
