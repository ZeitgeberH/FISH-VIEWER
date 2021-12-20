# -*- coding: utf-8 -*-
import pandas as pd
import os
import sys
import json
import tqdm
# import requests  ## if we query xml
# import xmltodict
from gr_allenFISH_utils import (
    check_gene_cached,
    check_gene_cached_images,
    check_cache_AtlasImage,
    # check_cache_AtlasDataFrame,
    load_cached_gene,
    load_cached_gene_images,
    load_cached_atlas_image,
    download_and_cache,
    download_and_cache_image,
    download_and_cache_AtlasImage,
    download_and_cache_AtlasImage_Annotated,
    Download_cache_AtlasDataFrame,
    extractStructureSVG,
)

# from brainrender import base_dir
from gr_io import request, fail_on_no_connection
# from brainrender.actors import Volume

from gr_sys_utils import base_dir

import pdb

class GeneRenderAPI:
    voxel_size = 200  # um
    grid_size = [58, 41, 67]  # number of voxels along each direction
    sectionPlaneDict = {1:'coronal',2:'sagittal'}
    DOWNSAMPLE_FACTOR = 4
    QUALITY_FACTOR = 50
    
    aba_mouseBrain_atlas = 'http://api.brain-map.org/api/v2/structure_graph_download/1.json'
    all_genes_url = (
        "http://api.brain-map.org/api/v2/data/query.json?criteria="
        + "model::Gene,"
        + "rma::criteria,products[abbreviation$eq'DevMouse'],"
        + "rma::options,[tabular$eq'genes.id','genes.acronym+as+gene_symbol','genes.name+as+gene_name',"
        + "'genes.entrez_id+as+entrez_gene_id','genes.homologene_id+as+homologene_group_id'],"
        + "[order$eq'genes.acronym']"
        + "&num_rows=all&start_row=0"
    )

    gene_experiments_url = (
        "http://api.brain-map.org/api/v2/data/query.json?criteria=model::SectionDataSet,"
        + "rma::criteria,[failed$eq'false'],products[abbreviation$eq'Mouse'],genes[acronym$eq-GENE_SYMBOL-]"
    )
    
    ## gridded data
    download_url = "http://api.brain-map.org/grid_data/download/EXP_ID?include=energy,intensity,density"
    
    ## xml query. more complex than json
    # download_url_full = "http://api.brain-map.org/api/v2/data/query.xml?criteria=model::SectionImage,rma::criteria,[data_set_id$eq-GENE_EXP_ID-]"
    
    ## section list
    download_url_full = "http://api.brain-map.org/api/v2/data/query.json?num_rows=all&criteria=model::SectionImage,rma::criteria,[data_set_id$eq-GENE_EXP_ID-]"
  
    ## thumbnail settings for nissl
    download_url_imageNissl_D = "http://api.brain-map.org/api/v2/image_download/-IMAGE_ID-?downsample=-DOWNSAMPLE_FACTOR-&quality=-QUALITY_FACTOR-"
    
    ## thumbnail settings for expression
    download_url_imageEXP_D = "http://api.brain-map.org/api/v2/image_download/-IMAGE_ID-?downsample=-DOWNSAMPLE_FACTOR-&quality=-QUALITY_FACTOR-&view=expression"
    
    ## section imagedata of Nissl-full-image
    download_url_imageNissl_F = "http://api.brain-map.org/api/v2/image_download/-IMAGE_ID-"
    ## section imagedata of EXPRESSION-full-image
    download_url_imageEXP_F = "http://api.brain-map.org/api/v2/image_download/-IMAGE_ID-?&view=expression"

    ##For each Structure in the input list, locate the closest image and (x,y) location in SectionDataSet 68545324:
    closed_image = "http://api.brain-map.org/api/v2/structure_to_image/-EXP_ID-.json?structure_ids=-STRUCTURE_ID-"

    ### Image-to_Atlas
    ##For Atlas 1 (P56 mouse), find the closest annotated SectionImage and (x,y) location as defined by a seed SectionImage and seed (x,y) location.
    closed_atlas1Image = "http://api.brain-map.org/api/v2/image_to_atlas/-IMAGE-ID.json?x=-X-&y=-Y-&atlas_id=1"
    
    gene_expression_cache = base_dir / "GriddedGeneExpressionCache"
    fish_images_cache = base_dir / "FishImagesCache"
    mouseAtlas_cache = base_dir / "ABA_mouseAtlas"
    adultMouse_atlas_path = mouseAtlas_cache /'aba_mouseAtlas.json' 
    gene_name = None
    gene_exp_id = None
    atlas_df = []
    
    def __init__(self):
        # Get metadata about all available genes
        self.genes = None  # when necessary gene data can be downloaded with self.get_all_genes
        self.gene_expression_cache.mkdir(exist_ok=True)
        self.fish_images_cache.mkdir(exist_ok=True)
        self.mouseAtlas_cache.mkdir(exist_ok=True)
        self.load_atlas_df()

    def getAtlasImagesAnnoationList(self, atlasID, GraphicGroupLabel_id=28):
        '''
        Returns list of atlasImage list with annotation
        -------
        data : TYPE: list of str (atlas image numbers)
            Atlases that have AtlasImages annotated with Structure boundaries, and the relevant GraphicGroupLabel.ids
            Adult Mouse, 3D Coronal:Atlas ID 602630314, GraphicGrouLabels: [28]
            Mouse P56, Coronal: Atlas ID 1, GraphicGrouLabels: [28,159226751]
            Mouse P56, Sagitall: Atlas ID 2, GraphicGrouLabels: [28,159226751]
        '''
        ##### all availabe list
        #url = "http://api.brain-map.org/api/v2/data/query.json?criteria=model::Atlas,rma::include,graphic_group_labels[name$il'Atlas*'],rma::options[only$eq'atlases.id,atlases.name,graphic_group_labels.id']"
        url = "http://api.brain-map.org/api/v2/data/query.csv?criteria=model::AtlasImage,rma::criteria,atlas_data_set(atlases[id$eq-ATLAS_ID-]),graphic_objects(graphic_group_label[id$eq-GRIPHIC_GROUP_LABEL_ID-]),rma::options[tabular$eq'sub_images.id'][order$eq'sub_images.id']&num_rows=all&start_row=0"
        url = url.replace("-ATLAS_ID-",str(atlasID)).replace("-GRIPHIC_GROUP_LABEL_ID-",str(GraphicGroupLabel_id))
        atlasImageList = request(url).content.decode().split()[1:]
        atlasImageList =[int(j) for j in atlasImageList]
        return atlasImageList,atlasID, GraphicGroupLabel_id

    def getAtlasImageAnnotated(self, atlasImageID, atlasID,grahpicGroupID=28):
        '''
        Download and caching annotaed atlas

        Parameters
        ----------
        atlasImageID : TYPE, int
            DESCRIPTION.
        atlasID : TYPE, int
            DESCRIPTION.
        grahpicGroupID : TYPE, int, optional
            DESCRIPTION. The default is 28.

        Returns
        -------
        None.

        '''
        url = "http://api.brain-map.org/api/v2/atlas_image_download/-ATLAS_IMG_ID-?annotation=true&atlas=-ATLAS-ID-"
        url = url.replace('-ATLAS_IMG_ID-',str(atlasImageID)).replace('-ATLAS-ID-',str(atlasID))
        download_and_cache_AtlasImage_Annotated(url, self.mouseAtlas_cache, atlasImageID, atlasID, grahpicGroupID)
        
    def getAtlasStructureDataFrame(self, atlasID,GraphicGroupLabel_id=28, use_cache=True):
        df_name = os.path.join(self.mouseAtlas_cache,'atlas_ID_'+str(atlasID)+'.pkl')
        if use_cache:
            if os.path.exists(df_name):
                cache=True
            else:
                cache=False
        else:
            cache=False
        if cache:
            print('Loading atlas dataframe from cache!')
            atlas_df = pd.read_pickle(df_name)
        else:
            atlasImageList,atlasID, GraphicGroupLabel_id = self.getAtlasImagesAnnoationList(atlasID,GraphicGroupLabel_id=28)
            url = "http://api.brain-map.org/api/v2/svg_download/-ALLAS_IMAGE_ID-?groups=-GRIPHIC_GROUP_LABEL_ID-"
            atlas_df = []
            for j in atlasImageList:
                print(f'image {j}')
                url = url.replace("-ALLAS_IMAGE_ID-",str(j)).replace("-GRIPHIC_GROUP_LABEL_ID-",str(GraphicGroupLabel_id))
                atlas_df.append(Download_cache_AtlasDataFrame(url, j, atlasID,GraphicGroupLabel_id))
            atlas_df = pd.concat(atlas_df)
            atlas_df.to_pickle(df_name)
        return atlas_df
        
    def getAtlasStructrueImage(self,strAcronym, atlasImageID, atlasID,GraphicGroupLabel_id=28):
        '''
        Extract a structure from atlas image

        Parameters
        ----------
        strAcronym : TYPE
            DESCRIPTION.
        atlasImageID : TYPE
            DESCRIPTION.
        atlasID : TYPE
            DESCRIPTION.
        GraphicGroupLabel_id : TYPE, optional
            DESCRIPTION. The default is 28.

        Returns
        -------
        list
            DESCRIPTION.

        '''
        url = "http://api.brain-map.org/api/v2/svg_download/-ALLAS_IMAGE_ID-?groups=-GRIPHIC_GROUP_LABEL_ID-"
        url = url.replace("-ALLAS_IMAGE_ID-",str(atlasImageID)).replace("-GRIPHIC_GROUP_LABEL_ID-",str(GraphicGroupLabel_id))
        q = self.get_structure_info(strAcronym)
        if len(q) > 0:
            strID = q['id'].values
        else:
            print(f'Structure {strAcronym} not found in atlas image {atlasImageID}')
            return []
        strSVG = extractStructureSVG(url, strID)
        if len(strSVG) >0:
            filePath = os.path.join(self.mouseAtlas_cache,strAcronym+'_'+str(atlasImageID)+'.svg')
            with open(filePath, 'w') as file:
                file.write(strSVG)
            
    def getAtlasBoundaryImage(self,atlasImageID, atlasID,GraphicGroupLabel_id=28):
        '''
        Download the structure boundary annotations (GraphicGroupLabel.id=28) for an AtlasImage (id=100960033) as a file (.svg):
        GraphicGroupLabels: https://help.brain-map.org/display/api/Atlas+Drawings+and+Ontologies
        http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Atlas,rma::include,graphic_group_labels[name$il%27Atlas*%27],rma::options[only$eq%27atlases.id,atlases.name,graphic_group_labels.id%27]
        Mouse P56, Coronal: Atlas ID 1, GraphicGrouLabels: [28,159226751]
        Mouse P56, Sagitall: Atlas ID 2, GraphicGrouLabels: [28,159226751]
        Parameters
        ----------
        imageID : TYPE
            DESCRIPTION.
        x : TYPE
            DESCRIPTION.
        y : TYPE
            DESCRIPTION.
        Returns
        -------
        None.

        '''
        fileExist = check_cache_AtlasImage(self.mouseAtlas_cache, atlasImageID, atlasID,GraphicGroupLabel_id)
        fileName = []
        if fileExist:
            fileName = load_cached_atlas_image(self.mouseAtlas_cache,atlasImageID, atlasID,GraphicGroupLabel_id)
        else:
            url = "http://api.brain-map.org/api/v2/svg_download/-ALLAS_IMAGE_ID-?groups=-GRIPHIC_GROUP_LABEL_ID-"
            url = url.replace("-ALLAS_IMAGE_ID-",str(atlasImageID)).replace("-GRIPHIC_GROUP_LABEL_ID-",str(GraphicGroupLabel_id))
            download_and_cache_AtlasImage(url, self.mouseAtlas_cache, atlasImageID, atlasID,GraphicGroupLabel_id)
            fileName = load_cached_atlas_image(self.mouseAtlas_cache,atlasImageID, atlasID,GraphicGroupLabel_id)
        return  fileName      
        
    def getExpTargetXY(self, expID, structureAcronym):
        '''
        Structure-To-Image
        For a target structure, find the closest SectionImage within an experiment and (x,y) location as defined by the centroid of the Structure.
        Parameters
        ----------
        expID : TYPE: Int
            Experiment ID.
        structureID : TYPE: Int
            Structure ID for a brain structure. Use function "get_structure_info()" to query with the structure's
            acyonmy as arugment.
        Returns
        -------
        Json data:
            section_image_id: Closest SectionImage to the 3-D centroid of the structure.
            section_number: Section number of the closest SectionImage.
            x: Closest x pixel coordinate of the closest SectionImage.
            y: Closest y pixel coordinate of the closest SectionImage.

        '''
        #closed_image = "http://api.brain-map.org/api/v2/structure_to_image/-EXP_ID-.json?structure_ids=-STRUCTURE_ID-"
        
        q = self.get_structure_info(structureAcronym)
        if len(q) > 0:
            structureID = q['id'].values
            print(f"{structureAcronym} ({(q['name'].values)[0]}) structure id: {structureID[0]}")
            url = self.closed_image.replace("-EXP_ID-", str(expID)).replace("-STRUCTURE_ID-", str(structureID[0]))
            data = request(url).json()["msg"][0]["image_sync"]
            return data
        else:
            print('data not found')
            return []
    
    def synImageToAtlas(self, atlasID, imageID, seedLocation):
        baseUrl = "http://api.brain-map.org/api/v2/image_to_atlas/-SECTION_IMAGE_ID-.json?x=-X-&y=-Y-&atlas_id=-ATLAS_ID-"
        try:
            url = baseUrl.replace("-ATLAS_ID-", str(atlasID)).replace("-X-", str(seedLocation[0])).replace("-Y-", str(seedLocation[1]))
            url = url.replace("-SECTION_IMAGE_ID-", str(imageID))
            return request(url).json()["msg"]["image_sync"]
        except:
            return []
        
        
    def load_atlas_df(self):
        if len(self.atlas_df) < 1:
            self.atlas_df = self.get_aba_mouseAtlas()
            
    def get_structure_info(self, structureName):
        '''        
        From aba atlas json file get Structure information
        Parameters
        ----------
        structureName : TYPE:str
            Allen Acyronym for a structure

        Returns
        -------
        q : TYPE: pandas dataframe
            Information about structure if exist
        '''
        df = self.atlas_df.copy()
        df['acronym']=[s.lower() for s in df['acronym']]
        q = df[df['acronym']==structureName.lower()]
        if len(q) >0:
            return q
        else:
            print(f'Structure {structureName} not found')
            return []
            
    @fail_on_no_connection
    def get_aba_mouseAtlas(self):
        """
        Load or download the "Mouse Brain Atlas" ontology as a hierarchically structured json file
        http://help.brain-map.org/display/api/Downloading+an+Ontology%27s+Structure+Graph
        """
        
        if os.path.exists(self.adultMouse_atlas_path):
            with open(self.adultMouse_atlas_path) as json_file:
                aba_json = json.load(json_file)
            # print('openning saved atlas')
        else:
            aba_json = request(self.aba_mouseBrain_atlas).json()['msg']
            with open(self.adultMouse_atlas_path, 'w') as outfile:
                json.dump(aba_json, outfile)
        root = [];
        out = [];
        aba_df = pd.DataFrame(self.flattenNestedJson(aba_json, root, out))
        return aba_df
    
    def flattenNestedJson(self, njs, root, out):
        """
        Flatten a nested json file

        Parameters
        ----------
        njs : TYPE: nest json
        root : TYPE: []
            DESCRIPTION.
        out : TYPE:[]
            output dictionary.

        Returns
        -------
        out : TYPE: list
            output dictionary.

        """
        for d in njs:
            if d['children'] == []:
                out.append(d)
            else:
                self.flattenNestedJson(d['children'], d, out)            
        return out
    
    @fail_on_no_connection
    def get_all_genes(self):
        """
        Download metadata about all the genes available in the Allen gene expression dataset
        """
        res = request(self.all_genes_url)
        return pd.DataFrame(res.json()["msg"])

    def get_gene_id_by_name(self, gene_name):
        self.gene_name = self.gene_name or gene_name
        if self.genes is None:
            self.genes = self.get_all_genes()

        if gene_name not in self.genes.gene_symbol.values:
            print(
                f"Gene name {gene_name} doesnt appear in the genes dataset, nothing to return\n"
                + "You can search for you gene here: https://mouse.brain-map.org/"
            )
            return None
        else:
            return int(
                self.genes.loc[self.genes.gene_symbol == gene_name].id.values[
                    0
                ]
            )

    def get_gene_symbol_by_id(self, gene_id):
        if self.genes is None:
            self.genes = self.get_all_genes()

        return self.genes.loc[
            self.genes.id == str(gene_id)
        ].gene_symbol.values[0]

    @fail_on_no_connection
    def get_gene_experiments(self, gene):
        """
        Given a gene_symbol it returns the list of ISH
        experiments for this gene as grided data

        :param gene_symbol: str
        """
        url = self.gene_experiments_url.replace("-GENE_SYMBOL-", gene)
        data = request(url).json()["msg"]

        if not len(data):
            print(f"No experiment found for gene {gene}")
            return None
        else:
            return [d["id"] for d in data]

    @fail_on_no_connection
    def get_gene_experiments2(self, gene):
        """
        Given a gene_symbol it returns the list of ISH
        experiments for this gene and full metaData

        :param gene_symbol: str
        """
        url = self.gene_experiments_url.replace("-GENE_SYMBOL-", gene)
        try:
            data = request(url).json()["msg"]
        except:
            print('API not responding for url: ', url)
            return
        if not len(data):
            print(f"No experiment found for gene {gene}")
            return None, None
        else:
            try: ## may fail and return internet error!
                return [d["id"] for d in data], data
            except:
                return None, None

    @fail_on_no_connection
    def get_gene_experiments_imageList(self, geneExpID):
        """
        Given a gene_experiment ID it returns the list of ISH
        imageIDs for this experiment

        :param geneExpID: str
        """
        url = self.download_url_full.replace("-GENE_EXP_ID-", str(geneExpID))

        # dict_data = xmltodict.parse(requests.get(url))  ## for xml query
        data = request(url).json()["msg"]
        if not len(data):
            print(f"No experiment found for gene {gene}")
            return None
        else:
            # nImages = int(dict_data['Response']['@total_rows'])
            # section_image = dict_data['Response']['section-images']['section-image'] ## list of OrderDict
            # return section_image            
            # nImages = len(data)
            imgdict = sorted({img['id']:img['section_number'] for img in data})
            imgIdx = {img['id']:idx for idx,img in enumerate(data)}
            sorted_data =[data[imgIdx[j]] for j in imgdict]
            section_number = [d['section_number'] for d in sorted_data]
            return imgdict, section_number, sorted_data

    def getBaseUrl(self,Expression=True, FullImage=False,D=4, Q=50):
        if FullImage:
            if Expression:
                baseUrl = self.download_url_imageEXP_F
            else:
                baseUrl = self.download_url_imageNissl_F                              
        else:
            if Expression:
                baseUrl = self.download_url_imageEXP_D
            else:
                baseUrl = self.download_url_imageNissl_D
            baseUrl = baseUrl.replace("-DOWNSAMPLE_FACTOR-", str(D))
            baseUrl = baseUrl.replace("-QUALITY_FACTOR-", str(Q))
        return baseUrl
    
    @fail_on_no_connection
    def download_gene_experiments_imageData(self, gene, expID, imageIDs,sectionNumbers, Expression=True, FullImage=False,D=4, Q=50):
        """
        Given an imageID it returns the jpeg image file

        :param imageIDs: list of numbers (image ID)
        :param section_number: list of numbers in the order of sectioning
        :param Expression: get Nissle or Expression value
        :param FullImage: get full image or downsampled and low-quality images (for quick browsing)
        """
        # download_url_imageEXP_D = "http://api.brain-map.org/api/v2/image_download/-IMAGE_ID-?downsample=-DOWNSAMPLE_FACTOR-&quality=-QUALITY_FACTOR-&view=expression"
       
        baseUrl = self.getBaseUrl(Expression, FullImage,D, Q)
         
        for j, sn in tqdm.tqdm(zip(imageIDs, sectionNumbers), total=len(imageIDs)):
            url = baseUrl.replace("-IMAGE_ID-", str(j))
            # print(f'Downlaoding image {j}')
            download_and_cache_image(url,os.path.join(self.fish_images_cache, f"{gene}-{expID}" ),sn,j,\
                                     Expression, FullImage,D,Q)
            
        
    @fail_on_no_connection
    def download_gene_data(self, gene):
        """
        Downloads a gene's data from the Allen Institute
        Gene Expression dataset and saves to cache.
        See: http://help.brain-map.org/display/api/Downloading+3-D+Expression+Grid+Data

        :param gene: int, the gene_id for the gene being downloaded.
        """
        # Get the gene's experiment id
        exp_ids = self.get_gene_experiments(gene)

        if exp_ids is None:
            return

        # download experiment data
        for eid in exp_ids:
            print(f"Downloading data for {gene} - experiment: {eid}")
            url = self.download_url.replace("EXP_ID", str(eid))
            download_and_cache(
                url, os.path.join(self.gene_expression_cache, f"{gene}-{eid}")
            )

    def get_gene_data(self, gene, exp_id, use_cache=True, metric="energy"):
        """
        Given a list of gene ids
        """
        self.gene_name = self.gene_name or gene

        # Check if gene-experiment cached
        if use_cache:
            cache = check_gene_cached(self.gene_expression_cache, gene, exp_id)
        else:
            cache = False

        if not cache:  # then download it
            self.download_gene_data(gene)
            cache = check_gene_cached(self.gene_expression_cache, gene, exp_id)
            if not cache:
                raise ValueError(  # pragma: no cover
                    "Something went wrong and data were not cached"
                )

        # Load from cache
        data = load_cached_gene(cache, metric, self.grid_size)

        if sys.platform == "darwin":
            data = data.T

        return data

    def griddata_to_volume(
        self,
        griddata,
        min_quantile=None,
        min_value=None,
        cmap="bwr",
    ):
        """
        Takes a 3d numpy array with volumetric gene expression
        and returns a vedo.Volume.isosurface actor.
        The isosurface needs a lower bound threshold, this can be
        either a user defined hard value (min_value) or the value
        corresponding to some percentile of the gene expression data.

        :param griddata: np.ndarray, 3d array with gene expression data
        :param min_quantile: float, percentile for threshold
        :param min_value: float, value for threshold
        """
        return Volume(
            griddata,
            min_quantile=min_quantile,
            voxel_size=self.voxel_size,
            min_value=min_value,
            cmap=cmap,
            name=self.gene_name,
            br_class="Gene Data",
        )
    
    @fail_on_no_connection    
    def get_gene_experiments_imageData(self, gene, exp_id, use_cache=True, Expression=False, FullImage=False,D=4, Q=50):
        """
        Given a gene and expID, return list of sections image data
        :param Expression: return expression data if true otherwise ISH data
        :param FullImage: return fullimage if true otherwise downsampling and lower-quatlity (using system default)
        """
        # Check if gene-experiment cached
        imageIDs, sectionNumbers, expMetaData = self.get_gene_experiments_imageList(exp_id) ## query imageIDs
        if use_cache:
            
            cache = check_gene_cached_images(self.fish_images_cache, gene, exp_id,Expression, FullImage, D, Q, len(imageIDs))
            # print(f'Availabe cache:{cache}')
        else:
            cache = False
            print(f'cache not availabe')

        if not cache:  # then download it

            if FullImage:
                print(f'Downloading images {imageIDs} with fullImage with Expression {Expression}')
            else:
                print(f'Downloading images:\n {imageIDs} \ndownsampling factor {self.DOWNSAMPLE_FACTOR} with quality {self.QUALITY_FACTOR}')                
            self.download_gene_experiments_imageData(gene, exp_id,imageIDs, sectionNumbers, Expression, FullImage,D, Q) ## download images
            cache = check_gene_cached_images(self.fish_images_cache, gene, exp_id,Expression, FullImage, D, Q,len(imageIDs))
            if not cache:
                raise ValueError(  # pragma: no cover
                    "Something went wrong and data were not cached"
                )

        # Load from cache
        data = load_cached_gene_images(cache,Expression)

        if sys.platform == "darwin":
            data = data.T

        return data
 
    
    @fail_on_no_connection    
    def download_gene_section_imageData(self, gene, expID, imageID, sectionNumber, Expression=False, FullImage=False, D=4, Q=50):
        """
        Given an imageID, download it from ABI
        :param Expression: return expression data if true otherwise ISH data
        :param FullImage: return fullimage if true otherwise downsampling and lower-quatlity (using system default)
        """             
        baseUrl = self.getBaseUrl(Expression, FullImage,D, Q)
         
        url = baseUrl.replace("-IMAGE_ID-", str(imageID))
        download_and_cache_image(url,os.path.join(self.fish_images_cache, f"{gene}-{expID}" ),sectionNumber,imageID,\
                                 Expression, FullImage,self.DOWNSAMPLE_FACTOR,self.QUALITY_FACTOR)