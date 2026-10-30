"""Faithful wrapper around the authors' BioMapAI (DNN.py) omics -> 12 scores -> disease model.
Hyper-parameters copied from codes/AI/module/Cross_validation.py: epochs=500, batch=64, Adam(5e-4), L2(0.008);
ScoreLayer: legacy Adam(0.01), L2(0.0008), 200 epochs; then WeightsAdjust (1 - train MSE per score).
The authors' predict() standardises the test set with its own mean/sd; we keep that behaviour."""
import os, sys, importlib.util
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import tensorflow as tf
tf.config.threading.set_intra_op_parallelism_threads(1)
tf.config.threading.set_inter_op_parallelism_threads(1)

_spec = importlib.util.spec_from_file_location('DNN', '/mnt/project-files/biomap/BioMapAI/codes/AI/module/DNN.py')
DNN = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(DNN)

# silence keras progress bars
for _m in ('fit', 'predict', 'evaluate'):
    _orig = getattr(tf.keras.Model, _m)
    def _wrap(self, *a, _o=_orig, **k):
        k.setdefault('verbose', 0); return _o(self, *a, **k)
    setattr(tf.keras.Model, _m, _wrap)


def fit_predict(Xtr, Str, ytr, Xte, epochs=500, seed=0):
    """Xtr: train features, Str: train 12 scores, ytr: train labels, Xte: test features.
    Returns (prob_test, predicted_scores_test)."""
    tf.keras.backend.clear_session()
    tf.random.set_seed(seed); np.random.seed(seed)
    osm = DNN.OmicScoreModel(epochs=epochs, optimizer=tf.keras.optimizers.Adam(0.0005), batch_size=64,
                             kernel_regularizer=tf.keras.regularizers.L2(0.008))
    model, _ = osm.train(np.asarray(Xtr, 'float32'), np.asarray(Str, 'float32'))
    sl = DNN.ScoreLayer(optimizer=tf.keras.optimizers.legacy.Adam(0.01),
                        kernel_regularizer=tf.keras.regularizers.l2(0.0008), epochs=200)
    score_layer = sl.build_model(np.asarray(Str, 'float32'), np.asarray(ytr, 'float32'))
    adj = DNN.WeightsAdjust(model, np.asarray(Xtr, 'float32'), np.asarray(Str, 'float32'), score_layer).adjust_score_weight()
    sy = DNN.ScoreYModel(model, adj)
    s_te = sy.get_score(np.asarray(Xte, 'float32'))
    p_te = np.asarray(adj.predict(s_te)).ravel()
    return p_te, s_te
