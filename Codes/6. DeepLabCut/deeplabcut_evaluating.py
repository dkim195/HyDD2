# deeplabcut_extract_frames.py
import deeplabcut

config_path = "/faststorage/project/Hypothalamus/deeplab/Dlx-LD-2024-07-21/config.yaml"

deeplabcut.evaluate_network(
    config_path,
    plotting=True,
)