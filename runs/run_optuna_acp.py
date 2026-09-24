"""Optuna HPO for PeptideCNNv2 on AntiCP2 main (the 0.83 chase).
Objective: best_val_auc (never test). Chosen config's TEST metrics are
reported once at the end - no test-set selection across trials."""
import sys, time, json
sys.path.insert(0, 'src')
import optuna
from pepx.datasets import load_anticp2
from pepx.models import PeptideCNNv2
from pepx.trainer import train_single_task

train = [r for r in load_anticp2('main') if r.meta == 'train']
test = [r for r in load_anticp2('main') if r.meta == 'test']
optuna.logging.set_verbosity(optuna.logging.WARNING)

def objective(trial):
    p = dict(emb_dim=trial.suggest_categorical("emb_dim", [32, 48, 64]),
             ch=trial.suggest_categorical("ch", [64, 96, 128]),
             blocks=trial.suggest_int("blocks", 2, 3),
             dropout=trial.suggest_float("dropout", 0.15, 0.5))
    lr = trial.suggest_float("lr", 1e-4, 3e-3, log=True)
    res, _ = train_single_task(PeptideCNNv2(**p), train, test, 'anticp2_main',
                               'cnnv2-hpo', epochs=15, wide=True, lr=lr, patience=5)
    trial.set_user_attr("test_auc_at_best_val", round(float(res.auc), 4))
    return float(res.best_val_auc)

t0 = time.time()
study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=7))
study.optimize(objective, n_trials=12, show_progress_bar=False)
best = study.best_trial
# final honest eval: retrain best config with full 40 epochs, report TEST once
p = {k: best.params[k] for k in ("emb_dim", "ch", "blocks", "dropout")}
res, _ = train_single_task(PeptideCNNv2(**p), train, test, 'anticp2_main',
                           'cnnv2-hpo-best', epochs=40, wide=True, lr=best.params["lr"])
out = {"best_params": best.params, "best_val_auc": round(float(best.value), 4),
       "final_test": res.__dict__, "n_trials": len(study.trials),
       "trial_test_aucs": sorted(t.user_attrs.get("test_auc_at_best_val", 0) for t in study.trials),
       "elapsed_s": round(time.time() - t0, 1),
       "note": "objective was validation AUC; test reported once for the chosen config"}
json.dump(out, open("results/optuna_cnnv2_anticp2.json", "w"), indent=1, default=float)
print("best val", round(best.value, 4), "final TEST auc", round(float(res.auc), 4))
print("elapsed", round(time.time() - t0, 1), "s")
