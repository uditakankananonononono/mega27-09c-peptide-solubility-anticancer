"""pepx HTTP API (FastAPI): POST /score with {"sequences": [...]} returns the
tri-objective probabilities from the shipped TriNet weights; GET /health."""
from __future__ import annotations
from pathlib import Path
from fastapi import FastAPI
from pydantic import BaseModel, Field

from .alphabet import InvalidSequenceError, validate_sequence
from .cli import _load_trinet, score_sequences

ROOT = Path(__file__).resolve().parents[2]
app = FastAPI(title="pepx", version="0.2.0")
_model, _mu, _sd = _load_trinet(ROOT / "results" / "trinet.pt")

class ScoreRequest(BaseModel):
    sequences: list[str] = Field(min_length=1, max_length=512)

@app.get("/health")
def health():
    return {"status": "ok", "model": "trinet-v0.1"}

@app.post("/score")
def score(req: ScoreRequest):
    seqs = []
    for s in req.sequences:
        try:
            seqs.append(validate_sequence(s))
        except InvalidSequenceError as e:
            return {"error": str(e)}
    scores = score_sequences(_model, seqs, _mu, _sd)
    return {"results": [
        {"sequence": s,
         "p_anticancer": round(float(scores["acp"][i]), 4),
         "p_soluble": round(float(scores["sol"][i]), 4),
         "p_aggregation": round(float(scores["agg"][i]), 4)}
        for i, s in enumerate(seqs)]}
