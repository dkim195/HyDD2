# deeplabcut_extract_frames.py
import deeplabcut

config_path = "/faststorage/project/Hypothalamus/deeplab/Dlx-LD-2024-07-21/config.yaml"

deeplabcut.create_multianimaltraining_dataset(
    config_path,
    num_shuffles=1,
    net_type="dlcrnet_ms5",
)

deeplabcut.train_network(
    config_path,
    saveiters=10000,
    maxiters=50000,
    allow_growth=True,
)
