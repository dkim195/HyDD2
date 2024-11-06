# deeplabcut_extract_frames.py
import deeplabcut

project_name = "Dlx"
experimenter = "LD"
video_path = "/faststorage/project/Hypothalamus/Video"
config_path = deeplabcut.create_new_project(
    project_name,
    experimenter,
    [video_path],
    multianimal=True,
    copy_videos=True,
)

#config_path = "/faststorage/project/Hypothalamus/deeplab/Dlx-TK-2024-07-18/config.yaml"

# Simulate the yes input
import builtins
builtins.input = lambda _: 'yes'

deeplabcut.extract_frames(
    config_path,
    mode="automatic",
    algo="kmeans",
    userfeedback=False,
)

