# FISH-VIEWER
## GUI to download, visualize FISH images from Allen Brain intitute
This is a packaged (python) program for downloading, browsing and analyzing Allen brain institute’s FISH image data using their APIs (https://help.brain-map.org/display/api/Allen+Brain+Atlas+API) . It also provides information for queried genes as well as Allen brain atlas ontology as a reference.

## Quick binary release
Download the binary release to your PC, unzip the file into a folder anywhere. Find the exe file name “FISH_VIEW.exe” and double click it. You can set up a short-cut to your windows desktop by right click that file and select “send/ Desktop (create shortcut)”. The first time you run the file, it may possibly trigger an alert from your anti-virus software. You can safely ignore it. Afterwards, It will bring up the main GUI.

## Main functionality
### Downloading files

There are two modes for downloading images: single gene query or batch downloading. 
Single gene query is best for exploring gene of interest. Enter the gene symbol (not case-sensitive) at top right corner (behind the label “Gene”, under the column of “Value”) and press “enter” key. The middle panel would list available experiments for this gene. You can restrict the image planes by choose “View options / Plane / Sagittal” (or Coronal). For Loading image, double click the experiment. This would automatically download all images for this experiment. By default, it will download the full resolution images of FISH data and its associated expression mask. Images will be rendered once it is done downloading. If you specified a brain region (in the parameter tab) to synchronize, the image that contain that brain region will be automatically choosing. 
Batch downloading. If you have a list of genes, you can write down their gene symbol (one per line) in a plain text file and store it anywhere in your drive. Then Click “Tools/Batch downloading” in the menu would create a dialog box. In the dialog box, click “locate file” and use the file dialog box to locate the file you wrote. Choose the resolution and planes you are interested in, click OK button.  After a few seconds, a dialog box should pop-out with information (Click ‘show Details’) about experiments associated with the gene listed in your file. Click ‘OK’ to download these experiments; Click ‘Abort’ to cancel downloading. If you have a working internet connection, this will download these images into your computer. Full resolution images are quite large (h>7k , w>10k ), so it would take some time.
Downloaded images are cached in a folder named “geneRender” under your user folder, for example: “C:\Users\MH\geneRender\FishImagesCache”. 

### Browsing images

The Tab “Section images” is for browsing images. You can use mouse or keyboard (up/down arrow) to go over images.  You can always ‘Syn’ back to the brain region you are interested by clicking the “Syn” button in the parameter tab above.

The two image panels are linked in X and Y axis. Moving either image will syn the other. Left-Clicking anywhere in the image to pan the image. Right clicking or tuning middle wheel to zoom.
To select an image region of interest, Changing the mouse mode to “1 button” (right clicking image to bring-up the context menu; Mouse mode / 1 button) then drag a rectangle around the region of interest. Click “view all” to bring back to full image. 
The context menu also has “export” option to export the image ROI to various format. 

Double clicking the ISH image would overlay atlas image in the bottom panel. User can adjust the position of the atlas image with UIs in the parameter panel.

### Helpful information.

The bottom has two tabs. First tab is about genes you recently queried.  Each gene has a http link (under column ‘genecard url’ ) to its GeneCard webpage. Double clicking the link would open a browser page to show the webpage. The second tab has ABA ontology. Clicking the column header to sort.

### Image analyzing	

Defining expressing pixels: “…the expressing pixels are colorized based on the luminosity of the ISH image (0.21 R + 0.72 G + 0.07 B) over some small local neighborhood”
See discussions in Allen Brain Institute’s forum on this issue: https://community.brain-map.org/t/quantification-of-gene-expression-by-in-situ-hybridization-finding-and-using-the-raw-values/153/10

And a paper about quantification of ISH expression:
Mackenzie Englund et. al 2021. Comparing cortex-wide gene expression patterns between species in a common reference frame.  https://www.biorxiv.org/content/10.1101/2021.07.28.454203v1.full

See this notebook showing mapping color mask back to intensity:
https://github.com/ZeitgeberH/Expression-color-mask-to-expresssion-intensity

## Screenshot

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
