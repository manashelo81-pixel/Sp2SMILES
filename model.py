"""Spectrum -> SELFIES transformer with an auxiliary fingerprint head."""
import math

import numpy as np
import selfies as sf
import torch
import torch.nn as nn
from rdkit import Chem

PAD, BOS, EOS, UNK = 0, 1, 2, 3


class MzEmbed(nn.Module):
    """Sinusoidal features over a wide range of wavelengths (0.01 .. 2000 Da)
    so the network can resolve fractions of a Da, plus intensity."""
    def __init__(self, d_model, n_freq=64, min_wl=0.01, max_wl=2000.0):
        super().__init__()
        wl = torch.exp(torch.linspace(math.log(min_wl), math.log(max_wl), n_freq))
        self.register_buffer("freq", 2 * math.pi / wl)
        self.mlp = nn.Sequential(nn.Linear(2 * n_freq + 1, d_model), nn.ReLU(),
                                 nn.Linear(d_model, d_model))

    def forward(self, mz, inten):
        ang = mz.unsqueeze(-1) * self.freq
        return self.mlp(torch.cat([ang.sin(), ang.cos(), inten.unsqueeze(-1)], -1))


class SpecNet(nn.Module):
    def __init__(self, vocab_size, n_adducts, n_instruments=1, d_model=256, nhead=8, n_enc=4, n_dec=4,
                 ff=1024, dropout=0.1, fp_bits=2048, max_len=160):
        super().__init__()
        self.mz = MzEmbed(d_model)
        self.adduct_emb = nn.Embedding(n_adducts, d_model)
        self.inst_emb = nn.Embedding(n_instruments, d_model)
        self.ce_proj = nn.Linear(2, d_model)  # [collision energy / 100, missing flag]
        enc = nn.TransformerEncoderLayer(d_model, nhead, ff, dropout, batch_first=True, norm_first=True)
        self.enc = nn.TransformerEncoder(enc, n_enc, enable_nested_tensor=False)
        self.tok = nn.Embedding(vocab_size, d_model, padding_idx=PAD)
        self.pos = nn.Embedding(max_len, d_model)
        dec = nn.TransformerDecoderLayer(d_model, nhead, ff, dropout, batch_first=True, norm_first=True)
        self.dec = nn.TransformerDecoder(dec, n_dec)
        self.out = nn.Linear(d_model, vocab_size)
        self.fp = nn.Linear(d_model, fp_bits)

    def encode(self, mz, inten, pad, prec, adduct, inst, ce):
        peaks = self.mz(mz, inten)
        one = torch.ones_like(prec).unsqueeze(1)
        p0 = (self.mz(prec.unsqueeze(1), one) + self.adduct_emb(adduct).unsqueeze(1)
              + self.inst_emb(inst).unsqueeze(1) + self.ce_proj(ce).unsqueeze(1))
        x = torch.cat([p0, peaks], 1)
        key_pad = torch.cat([torch.zeros_like(pad[:, :1]), pad], 1)
        mem = self.enc(x, src_key_padding_mask=key_pad)
        keep = (~key_pad).unsqueeze(-1).float()
        pooled = (mem * keep).sum(1) / keep.sum(1).clamp(min=1)
        return mem, key_pad, self.fp(pooled)

    def decode(self, tokens, mem, key_pad):
        L = tokens.size(1)
        y = self.tok(tokens) + self.pos(torch.arange(L, device=tokens.device))
        causal = torch.triu(torch.ones(L, L, dtype=torch.bool, device=tokens.device), 1)
        h = self.dec(y, mem, tgt_mask=causal, tgt_key_padding_mask=(tokens == PAD),
                     memory_key_padding_mask=key_pad)
        return self.out(h)


def make_batch(items, device, peak_dropout=0.0, rng=None):
    """items: list of (peaks float32 (n,2) with intensity 0..1, precursor_mz, adduct_id,
    instrument_id, collision_energy or None)."""
    rng = rng or np.random
    n = max(len(x[0]) for x in items)
    B = len(items)
    mz = np.zeros((B, n), np.float32)
    it = np.zeros((B, n), np.float32)
    pad = np.ones((B, n), bool)
    for i, (p, *_) in enumerate(items):
        if peak_dropout > 0 and len(p) > 5:
            keep = rng.random(len(p)) > peak_dropout
            if keep.sum() >= 3:
                p = p[keep]
        k = len(p)
        mz[i, :k] = p[:, 0]
        it[i, :k] = np.sqrt(p[:, 1])
        pad[i, :k] = False
    prec = np.array([x[1] for x in items], np.float32)
    add = np.array([x[2] for x in items], np.int64)
    inst = np.array([x[3] for x in items], np.int64)
    ce = np.array([[0.0, 1.0] if x[4] is None else [min(x[4], 200.0) / 100.0, 0.0] for x in items], np.float32)
    t = lambda a: torch.from_numpy(a).to(device)
    return t(mz), t(it), t(pad), t(prec), t(add), t(inst), t(ce)


@torch.no_grad()
def generate(model, mem, key_pad, max_len=130, temperature=0.0):
    B = mem.size(0)
    toks = torch.full((B, 1), BOS, dtype=torch.long, device=mem.device)
    done = torch.zeros(B, dtype=torch.bool, device=mem.device)
    for _ in range(max_len):
        logits = model.decode(toks, mem, key_pad)[:, -1].float()
        logits[:, [PAD, BOS, UNK]] = -1e9
        if temperature <= 0:
            nxt = logits.argmax(-1)
        else:
            nxt = torch.multinomial(torch.softmax(logits / temperature, -1), 1).squeeze(1)
        nxt = torch.where(done, torch.full_like(nxt, PAD), nxt)
        toks = torch.cat([toks, nxt.unsqueeze(1)], 1)
        done |= nxt == EOS
        if done.all():
            break
    return toks[:, 1:].cpu().numpy()


def tokens_to_mol(ids, itos):
    """Token ids -> (flat canonical SMILES, mol) or (None, None) if invalid."""
    toks = []
    for i in ids:
        if i == EOS or i == PAD:
            break
        toks.append(itos[int(i)])
    if not toks:
        return None, None
    try:
        smi = sf.decoder("".join(toks))
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None, None
        flat = Chem.MolToSmiles(mol, isomericSmiles=False)
        return flat, Chem.MolFromSmiles(flat)
    except Exception:
        return None, None


def load_model(path, device):
    ck = torch.load(path, map_location=device, weights_only=False)
    model = SpecNet(**ck["model_args"]).to(device)
    model.load_state_dict(ck["state_dict"])
    model.eval()
    return model, ck
