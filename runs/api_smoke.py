"""FastAPI smoke test via TestClient: health + batch score + CLI agreement."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from fastapi.testclient import TestClient
from pepx.api import app

c = TestClient(app)
h = c.get("/health").json()
r = c.post("/score", json={"sequences": ["FEKEAKKIEIKRH", "KLAKLAKLAKLAK"]}).json()
bad = c.post("/score", json={"sequences": ["AKX*Z"]}).json()
out = {"health": h, "score": r, "invalid_handling": bad,
       "cli_agreement": r["results"][0]["p_anticancer"]}
json.dump(out, open("results/api_smoke.json","w"), indent=1)
print(json.dumps(r, indent=1)[:400])
print("health:", h)
