import os
import numpy as np
import zipfile
import io
import sys
import shutil
from PIL import Image
from gr_sys_utils import get_subdirs, listdir
from gr_io import request, check_file_exists
import pdb
import glob
import pandas as pd
from biothings_client import get_client
import scipy.cluster.vq as scv
geclient = get_client('gene')
import ncbi.datasets
import matplotlib.cm as cm
# Start a Datasets gene API instance
ncbi_api_client = ncbi.datasets.ApiClient()
ncbi_gene_instance = ncbi.datasets.GeneApi(ncbi_api_client)


## Expression Energy

#The expression mask image display highlights those cells that have the highest probability of gene expression using a heat map color scale (from low/blue to high/red).
#https://help.brain-map.org/display/mousebrain/In+Situ+Hybridization+%28ISH%29+Data

## Quantification of gene expression by in situ hybridization: finding and using the raw values
##https://community.brain-map.org/t/quantification-of-gene-expression-by-in-situ-hybridization-finding-and-using-the-raw-values/153

'''
Once the segmentation algorithm has identified expressing vs non-expressing pixels, the expressing pixels are colorized based on the luminosity of the ISH image (0.21 R + 0.72 G + 0.07 B) over some small local neighborhood. On the web application it is then colorized using the jet colormap. As such this not a “quantity” of expression as you put it. Rather, we wanted to provide a view of the color intensity in the expression mask.

Segmentation takes into account intensity of the local background which may contributes to the overlap that you are observing.

Lydia


'''
def colormap2arr(arr,cmap, cutoff=1.0):    
    # Modified from http://stackoverflow.com/questions/3720840/how-to-reverse-color-map-image-to-scalar-values/3722674#3722674
    N = 64
    gradient = np.zeros((N,4))
    gradient[1:,:]=cmap(np.linspace(0.0,1.0,N-1))
    if cmap==cm.jet: ## for jet color map, the last few entries are closer to origin
        maxRed=np.argwhere(gradient[:,0][::-1]>=cutoff)[0][0] 
        gradient = gradient[:-maxRed,:] ## remove those before mapping   
    # Reshape arr to N by 4
    arr2=arr.reshape((arr.shape[0]*arr.shape[1],arr.shape[2]))
    # Use vector quantization to shift the values in arr2 to the nearest point in
    # the code book (gradient).
    code,dist=scv.vq(arr2,gradient)
    # Scale the values so they are from 0 to 1.
    values=code.astype('float')/gradient.shape[0]
    # Reshape values back to original shape
    values=values.reshape(arr.shape[0],arr.shape[1])
    return values

def prepArr(arr):
    arr4c = np.ones((arr.shape[0],arr.shape[1],4))
    arr_copy = arr.copy()/255.0 ## normalize to range(0,1.0)
    arr4c[:,:,:3] = arr_copy ## return 4 channel array
    return arr4c

def rgb2intensity(arr,cmap=cm.jet):
    arr4c = prepArr(arr)    
    return colormap2arr(arr4c,cmap)

def check_gene_cached(cache_folder, gene_id, exp_id):
    """
    A gene is saved in a folder in cache_folder
    with gene_id-exp_id as name. If the folder doesn't
    exist the gene is not cached.

    :param cache_folder: str, path to general cache folder for all data
    :param gene_id: str name of gene
    :param exp_id: id of experiment
    """
    cache = [
        sub
        for sub in get_subdirs(cache_folder)
        if f"{gene_id}-{exp_id}" == os.path.basename(sub)
    ]
    if not cache:
        return False
    elif len(cache) > 1:
        raise ValueError("Found too many folders")
    else:
        return cache[0]

def check_gene_cached_images(cache_folder, gene_id, exp_id,Expression=False, FullImage=False, D=4, Q=50, nimages=1):
    """
    A gene is saved in a folder in cache_folder
    with gene_id-exp_id as name. If the folder doesn't
    exist the gene is not cached.

    :param cache_folder: str, path to general cache folder for all data
    :param gene_id: str name of gene
    :param exp_id: id of experiment
    """
    cache = [
        sub
        for sub in get_subdirs(cache_folder)
        if f"{gene_id}-{exp_id}" == os.path.basename(sub)
    ]
    if not cache:
        return False
    elif len(cache) > 1:
        raise ValueError("Found too many folders")
    else:

        
        cache_sub = [
            sub.split('\\')[-1]
            for sub in get_subdirs(cache[0])           
        ]
        # print(f'Found cache for {gene_id}-{exp_id} with subfolders: {cache_sub}')
        if FullImage:
            if 'F' in cache_sub:
                if Expression:
                    nImages = len(glob.glob(os.path.join(cache[0],'F/*_E.jpg')))

                else:
                    nImages = len(glob.glob(os.path.join(cache[0],'F/*_H.jpg')))
                if nImages > nimages-1:
                    return os.path.join(cache[0],'F')
                else:
                    return False
                    
            else:
                return False
        else:
            if 'DQ' in cache_sub:
                if Expression:
                    nImages = len(glob.glob(os.path.join(cache[0],'DQ/*_E.jpg')))

                else:
                    nImages = len(glob.glob(os.path.join(cache[0],'DQ/*_H.jpg')))
                if nImages > nimages-1:
                    return os.path.join(cache[0],'DQ')
                else:
                    return False
                    
            else:
                return False
            
def download_and_cache(url, cachedir):
    """
    Given a url to download a gene's ISH experiment data,
    this function download and unzips the data

    :param url: str, utl to download data
    :param cachedir: str, path to folder where data will be downloaded
    """
    # Get data
    req = request(url)

    # Create cache dir
    if not os.path.isdir(cachedir):
        os.mkdir(cachedir)

    # Unzip to cache dir
    z = zipfile.ZipFile(io.BytesIO(req.content))
    z.extractall(cachedir)

def download_and_cache_image(url, cachedir,index, imageID,Expression=True, FullImage=False,D=4,Q=50):
    """
    Given a url to download an image for gene's ISH experiment data with specificed quality
    :param url: str, utl to download data
    :param cachedir: str, path to folder where data will be downloaded
    :param D: downsampling factor
    :param Q: quality factor
    """

    # Create cache dir
    if not os.path.isdir(cachedir):
        os.mkdir(cachedir)
    if FullImage:
        fileDir = os.path.join(cachedir,'F')
        if not os.path.isdir(fileDir):
            os.mkdir(fileDir)
    else:
        fileDir = os.path.join(cachedir,'DQ')
        if not os.path.isdir(fileDir):
            os.mkdir(fileDir)
                        
    #reqest image
    req = request(url)
    z = Image.open(io.BytesIO(req.content)) ## decode image (Jpeg Image file.Mode RGB)
    if Expression:
        fileName = os.path.join(fileDir,str(index)+'_'+str(imageID)+'_E.jpg') ## expression images
    else:
        fileName = os.path.join(fileDir,str(index)+'_'+str(imageID)+'_H.jpg') ## ISH images
    z.save(fileName)

def check_cache_AtlasImage(cachedir,AtlasImageID, atlasID,grahpicGroupID):
    """
    Given a url to download an image for an atlas

    """
    # Create cache dir
    if not os.path.isdir(cachedir):
        return False    
    groupDir = os.path.join(cachedir,str(grahpicGroupID))
    if not os.path.isdir(groupDir):
        return False
        
    atlasDir = os.path.join(groupDir,str(atlasID))
    if not os.path.isdir(atlasDir):
        return False                          
    fileName = os.path.join(atlasDir,str(AtlasImageID)+'.svg')
    return os.path.exists(fileName)

def download_and_cache_AtlasImage_Annotated(url, cachedir,AtlasImageID, atlasID, grahpicGroupID):
    if not os.path.isdir(cachedir):
        os.mkdir(cachedir)
    
    groupDir = os.path.join(cachedir,str(grahpicGroupID))
    if not os.path.isdir(groupDir):
        os.mkdir(groupDir)
        
    atlasDir = os.path.join(groupDir,str(atlasID))
    if not os.path.isdir(atlasDir):
        os.mkdir(atlasDir)
    #reqest image
    req = request(url)
    z = Image.open(io.BytesIO(req.content)) ## decode image (Jpeg Image file.Mode RGB)
    fileName = os.path.join(atlasDir,str(AtlasImageID)+'_annotation.jpg') ## ISH images
    z.save(fileName)
    
    
def download_and_cache_AtlasImage(url, cachedir,AtlasImageID, atlasID, grahpicGroupID):
    """
    Given a url to download an image for an atlas

    """
    # Create cache dir
    if not os.path.isdir(cachedir):
        os.mkdir(cachedir)
    
    groupDir = os.path.join(cachedir,str(grahpicGroupID))
    if not os.path.isdir(groupDir):
        os.mkdir(groupDir)
        
    atlasDir = os.path.join(groupDir,str(atlasID))
    if not os.path.isdir(atlasDir):
        os.mkdir(atlasDir)                           
    #reqest image
    svgText = request(url).text
    fileName = os.path.join(atlasDir,str(AtlasImageID)+'.svg') ## expression images
    print('downloading atlas image: '+str(AtlasImageID))
    with open(fileName, 'w') as file:
        file.write(svgText)
    print('Done!')

def load_cached_gene(cache, metric, grid_size):
    """
    Loads a gene's data from cache
    """
    files = [
        f for f in listdir(cache) if metric in f and not f.endswith(".mhd")
    ]
    if not files:
        return None
    if len(files) > 1:
        raise NotImplementedError("Deal with more than one file found")
    else:
        return read_raw(files[0], grid_size)

def load_cached_gene_images(cache,Expression):
    if Expression:
        surfix = 'E'
    else:
        surfix = 'H'
    print(f'cache gene expression\n: {cache}')
    imagesNames = glob.glob(os.path.join(cache,'*_'+surfix+'.jpg'))
    img_section=np.argsort([int(x.split('\\')[-1].split('_')[0]) for x in imagesNames])   
    return [imagesNames[j] for j in img_section]

def load_cached_atlas_image(cachedir,atlasImageID, atlasID,grahpicGroupID):
    groupDir = os.path.join(cachedir,str(grahpicGroupID))     
    atlasDir = os.path.join(groupDir,str(atlasID))
    fileName = os.path.join(atlasDir,str(atlasImageID)+'.svg')
    return fileName

def Download_cache_AtlasDataFrame(url,atlasImgID, atlasID,GraphicGroupLabel_id):
    '''
    Make a dataframe Identify all structures in current atlas images

    Parameters
    ----------
    url : TYPE: str
        http address of SVG atlas downloading from ABA atlas boundary images

    Returns
    -------
    list of int: structure_id

    '''
    svgText = request(url).text
    strLst = []
    strIdx = svgText.find('structure_id="') 
    if strIdx==-1:
        print('no structure found')
        return []
    while strIdx!=-1:              
        strIdx2 = svgText[strIdx+13:].find(' d="') ## len(structure_id=") is 13
        strLst.append(svgText[strIdx+14:strIdx+12+strIdx2])
        svgText = svgText[strIdx+13+strIdx2:]
        strIdx = svgText.find('structure_id="') ## looking for next structure
    strLst = [int(j) for j in strLst]
    df = pd.DataFrame(strLst, columns=['id'])
    df['atlasImgID'] = atlasImgID
    df['atlasID'] = atlasID
    df['GraphicGroupLabel_id'] = GraphicGroupLabel_id
    return df
           
def extractStructureSVG(url,idx):
    '''    
    Extract SVG fragment correcponding to parituclar structure
    Parameters
    ----------
    url : TYPE: str
        http address of SVG atlas downloading from ABA atlas boundary images
    idx : TYPE: int
        Anatomical structure ID.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    svgText = request(url).text
    strIdx = svgText.find(f'structure_id="{str(idx)}"')
    if strIdx==-1:
        print(f'no structure with id {idx} found')
        return []
    structureStrings = ''
    head = svgText[:svgText.find('<path id')]
    tail = '</g></g></svg>'
    tailIdx = svgText.find(f'structure_id="{str(idx)}"')
    pathIdx1 = svgText[strIdx-60:strIdx].find('<path id')
    pathIdx1 = strIdx-60+pathIdx1
    pathIdx2 = svgText[strIdx:].find('<path id')
    pathIdx2 = strIdx+pathIdx2
    structureStrings=structureStrings+svgText[pathIdx1:pathIdx2]
    while pathIdx2 < tailIdx:
        strIdx = svgText[pathIdx2:].find(f'structure_id="{str(idx)}"')
        if strIdx==-1:
            print(f'no structure with id {idx} found')
            return head+structureStrings+tail
        pathIdx1 = svgText[strIdx-60:strIdx].find('<path id') ## 60 characters backwards are just right!
        pathIdx1 = strIdx-60+pathIdx1
        pathIdx2 = svgText[strIdx:].find('<path id')
        pathIdx2 = strIdx+pathIdx2 
        structureStrings=structureStrings+svgText[pathIdx1:pathIdx2]
    return head+structureStrings+tail
    
# --------------------------------- Open .raw -------------------------------- #
@check_file_exists
def read_raw(filepath, grid_size):
    """
    reads a .raw file with gene expression data
    downloaded from the Allen atlas and returns
    a numpy array with the correct grid_size.
    See as reference:
        http://help.brain-map.org/display/mousebrain/API#API-Expression3DGridsz

    :param filepath: str or Path object
    """
    filepath = str(filepath)

    # Read bytes
    with open(filepath, "rb") as test:
        content = test.read()

    # Create np array and return
    data = np.frombuffer(content, dtype="float32").reshape(grid_size)

    if sys.platform == "darwin":
        data = data.T  # TODO figure out why this is necessary on Mac OS?

    return data

# def readJpegImage(f):  ## need pvvips and libvips. may not work with pyinstaller
#     image = pyvips.Image.new_from_file(f, access="sequential") 
#     mem_img = image.write_to_memory() 
#     imgnp=np.frombuffer(mem_img, dtype=np.uint8)#.reshape(image.height, image.width, 3)  
#     return imgnp

def getDefaultBrowerPath():
    from winreg import HKEY_CLASSES_ROOT, HKEY_CURRENT_USER, OpenKey, QueryValueEx    
    try:
        with OpenKey(HKEY_CURRENT_USER, r'SOFTWARE\Microsoft\Windows\Shell\Associations\UrlAssociations\http\UserChoice') as regkey:
            # Get the user choice
            browser_choice = QueryValueEx(regkey, 'ProgId')[0]
    
        with OpenKey(HKEY_CLASSES_ROOT, r'{}\shell\open\command'.format(browser_choice)) as regkey:
            # Get the application the user's choice refers to in the application registrations
            browser_path_tuple = QueryValueEx(regkey, None)
    
            # This is a bit sketchy and assumes that the path will always be in double quotes
            browser_path = browser_path_tuple[0].split('"')[1]
            return browser_path
    except Exception:
        log.error('Failed to look up default browser in system registry. Using fallback value.')
        return -1
        
def queryGene(geneSymbol, species='mouse'):
    '''
    return gene information from mygene.com

    Parameters
    ----------
    geneSymbol : TYPE
        DESCRIPTION.
    species : TYPE, optional
        DESCRIPTION. The default is 'mouse'.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    try:
        df = geclient.query('symbol:'+geneSymbol,species=species.lower())
        if df:
            return df['hits'][0]
        else:
            return []
    except:
        return []
   
def queryGene_ncbi(geneSymbol, taxon='mouse'):
    '''
    return gene information from NCBI database using ncbi API

    Parameters
    ----------
    geneSymbol : TYPE, str
        DESCRIPTION.
    taxon : TYPE, optional
        DESCRIPTION. The default is 'mouse'.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    gene_metadata = ncbi_gene_instance.gene_metadata_by_tax_and_symbol(symbols=[geneSymbol], taxon=taxon)
    if gene_metadata.genes[0].gene!=None:
        dict = {}
        gene = gene_metadata.genes[0].gene
        dict['description'] = gene['description']
        dict['entrez_id'] = gene['gene_id']
        dict['type'] = gene['type']
        dict['chromosomes'] = '_'.join(gene['chromosomes'])
        dict['number of transcripts'] = len(gene['transcripts'])
        return dict
    else:
        return []