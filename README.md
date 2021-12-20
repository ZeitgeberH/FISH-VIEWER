# FISH-VIEWER
## GUI to download, visualize FISH images from Allen Brain intitute

### FISH image and expression mask are auto synced.

![expression mask](/figure/Gad1_colorMask.PNG)

### FISH image and inferred expression intensity

![expression mask](/figure/Gad1_intensity.PNG)

## To make excutable
Create conda enviorment by:
```
conda create --name fishview --file fishview_env_spec-file.txt
```

After activate your enviorment,
```
conda activate fishview
```

Using Pyinstaller to make package:
```
pyinstaller --clean .\FISHViewer_folder.spec
```
