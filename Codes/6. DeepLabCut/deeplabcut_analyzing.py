# deeplabcut_extract_frames.py
import deeplabcut

config_path = "/faststorage/project/Hypothalamus/script/deep/Dlx-LD-2024-07-22/config.yaml"

new_video_path = ['/faststorage/project/Hypothalamus/script/deep/Dlx-LD-2024-07-22/new_video.mp4']

deeplabcut.analyze_videos(
    config_path,
    [video],
    auto_track=True,
)

#deeplabcut.convert_detections2tracklets(
#    config_path,
#    [video],
#    track_method="ellipse",
#)

#deeplabcut.stitch_tracklets(
#    config_path,
#    [video],
#    track_method="ellipse",
#    min_length=5,
#)

deeplabcut.create_labeled_video(
    config_path,
    [video],
    color_by="individual",
    keypoints_only=False,
    trailpoints=10,
    draw_skeleton=False,
    track_method="ellipse",
)