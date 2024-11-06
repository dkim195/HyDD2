import deeplabcut

#Create a New Project - do with GUI
project_name = "Dlx"
experimenter = "LD"
video_path = "/faststorage/project/Hypothalamus/Video"
config_path = deeplabcut.create_new_project(
    project_name,
    experimenter,
    [video_path],
    multianimal=False,
    copy_videos=True,
)

#Configure the Project  - do with GUI
    # Head, Neck, Tail + Objective
#Add new vidoes if necessary

#extract - do with GUI
deeplabcut.extract_frames(config_path, mode='automatic', algo='kmeans', crop=False)

#label -  do with GUI - manual
deeplabcut.label_frames(config_path)

#create training  - do with GUI - fast
deeplabcut.create_training_dataset(config_path)

#training  - do with GUI
deeplabcut.train_network(config_path)

#evaluating  - do with GUI 
deeplabcut.evaluate_network(config_path, plotting=True)

#analyze - Python
config_path = r"C:\Users\au731006\Downloads\Dlx-LD-2024-07-25\config.yaml"
video_path = r'C:\Users\au731006\Downloads\Dlx-LD-2024-07-25\videos\IMG_1531_2.mov'
video_path1 = r'C:\Users\au731006\Downloads\Dlx-LD-2024-07-25\videos\IMG_1530_2.mov'
output_dir = r'C:\Users\au731006\Downloads\Dlx-LD-2024-07-25\videos'

deeplabcut.analyze_videos(config_path, video_path, videotype='mov', save_as_csv=True, destfolder=output_dir)
#or in GUI
    #Tick - Filter predictions, Plot trajectories, Show trajectory plots, Save results(s) as csv


#Track and refine - Python
deeplabcut.convert_detections2tracklets(config_path, video_path, track_method="ellipse") #fast
deeplabcut.stitch_tracklets(config_path, video_path, track_method="ellipse", min_length=5)
deeplabcut.filterpredictions(config_path, video_path, track_method="ellipse")
deeplabcut.plot_trajectories(config_path, video_path, filtered=True, track_method="ellipse")

#create - haven't done this
deeplabcut.create_labeled_video(
    config_path,
    video_path,
    color_by="individual",
    keypoints_only=False,
    trailpoints=10,
    draw_skeleton=True,
    track_method="ellipse",
    filtered=True
)

#if refinement is necessary
deeplabcut.extract_outlier_frames(config_path, video_path)
deeplabcut.refine_labels(config_path)
deeplabcut.merge_datasets(config_path)
deeplabcut.create_training_dataset(config_path)
deeplabcut.train_network(config_path, shuffle=2)  # Using a new shuffle







