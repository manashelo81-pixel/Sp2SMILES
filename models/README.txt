Put your trained model and cleaned library here (or point the app to them elsewhere):

  model.pt      from 02_train.py   (run folder, e.g. run2/model.pt)
  dataset.pkl   from 01_prepare_data.py   (processed_all/dataset.pkl)

On the first start the app builds a fingerprint cache next to dataset.pkl (dataset.pkl.pool.pkl, a few minutes,
one time only).
