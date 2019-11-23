import os
import numpy as np
import deepchem as dc
from deepchem.molnet import load_pdbbind_grid
pdbbind_tasks, pdbbind_datasets, transformers = load_pdbbind_grid(split='random', subset='full')
train_data, valid_data, test_data = pdbbind_datasets
metric = dc.metrics.Metric(dc.metrics.pearson_r2_score)
model = dc.models.MultitaskRegressor(len(pdbbind_tasks), train_data.X.shape[1],dropouts=[.25],learning_rate=0.0003,weight_init_stddevs=[.1],batch_size=64, model_dir='pdbbind_test')

model.fit(train_data, nb_epoch=100)
train_score = model.evaluate(train_data, [metric], transformers)
valid_score = model.evaluate(valid_data, [metric], transformers)

feature = dc.feat.RdkitGridFeaturizer(voxel_width=16.0,feature_types=["ecfp", "splif", "hbond", "salt_bridge"],ecfp_power=9,splif_power=9,flatten=True)

grid = feature.featurize_complexes(['ligand2.sdf'], ['3lpt.pdb'])
print(grid)

# reload = dc.models.MultitaskRegressor(len(pdbbind_tasks), train_data.X.shape[1],dropouts=[.25],learning_rate=0.0003,weight_init_stddevs=[.1],batch_size=64, model_dir='pdbbind_test')
# print(reload.get_checkpoints())
# reload.load_from_pretrained(reload, model_dir='pdbbind_test')
# reload.evaluate(train_data, [metric], transformers)
# reload.evaluate(valid_data, [metric], transformers)