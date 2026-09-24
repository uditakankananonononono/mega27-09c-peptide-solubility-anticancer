import sys, time, json
sys.path.insert(0, 'src')
t0 = time.time()
from pepx.datasets import augmented_acp_train, load_anticp2
from pepx.models import PeptideCNNv2
from pepx.trainer import train_single_task
train = augmented_acp_train('main')
test = [r for r in load_anticp2('main') if r.meta == 'test']
res, prob = train_single_task(PeptideCNNv2(), train, test, 'anticp2_main',
                              'PeptideCNNv2+aug', epochs=40, wide=True)
print(res)
print('elapsed', round(time.time() - t0, 1), 's')
with open('results/cnnv2_aug_anticp2_main.json', 'w') as fh:
    json.dump(res.__dict__, fh, indent=2)
