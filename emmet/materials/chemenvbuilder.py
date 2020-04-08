import json
import os
from pprint import pprint
from pymatgen import Structure
import urllib as url
from pprint import pprint
from pydash import get

from pymongo import MongoClient
from tqdm import tqdm
# from tqdm.autonotebook import tqdm
from maggma.core import Builder
from maggma.core import Store
from maggma.builders.map_builder import MapBuilder
from maggma.stores import MongograntStore

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
import logging
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy, MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments

# Connecting to approx_neb database
approx_neb = MongograntStore("rw:mongodb03.nersc.gov/mp_jmlin_mv","approx_neb_jerrycopy")
approx_neb.connect()

# connecting to sandbox database
sandbox = MongograntStore("rw:mongodb03.nersc.gov/mp_jmlin_mv","sandbox")
sandbox.connect()

# renaming for convenience
source_collection = approx_neb
target_collection = sandbox

# changing key to "wf_uuid" since convention in approx_neb
source_collection.key = "wf_uuid"
target_collection.key = "wf_uuid"

class ChemEnvCoordEnvBuilder(MapBuilder):
    """
    given a doc (dict) from the approx_neb collection, analyzes the coordination environment of 
    the sites along path (start site, path sites, end site) and returns the result in a new doc (dict)
    using ChemEnv from pymatgen
    
    Args:
        source_collection[MongoStore]: collection specified by MongoStore() from maggma.stores, input collection
        target_collection[MongoStore]: collection specified by MongoStore() from maggma.stores, output collection
    """
    def analyze_struct_chemenv(self, s, elem):
        """helper function that given a structure (s) and a working ion (elem), analyzes the 
        coordination environment using ChemEnv and [WHAT TYPE OF STRATEGY] and return the output as a dict
        
        NOTE: currently using SimplestChemenvStrategy w/ distance_cutoff=1.4 and angle_cutoff=0.3"""
        # Setup the local geometry finder
        lgf = LocalGeometryFinder()
        # you can also save the logging to a file, just remove the comment
#         logging.basicConfig(#filename='chemenv_structure_environments.log',
#                             format='%(levelname)s:%(module)s:%(funcName)s:%(message)s',
#                             level=logging.DEBUG)
        lgf.setup_structure(structure=s)
        # Get the StructureEnvironments and only looking at working ion
        se = lgf.compute_structure_environments(maximum_distance_factor=1.41, only_atoms=[elem])
        strategy = SimplestChemenvStrategy(distance_cutoff=1.4, angle_cutoff=0.3)
        lse = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=se)
        analysis = {'structure':s.as_dict(),               # save results into output dictionary
                    'working_ion':elem,
                    'coord_env':lse.coordination_environments[0][0]['ce_symbol'],
                    'coord_num':lse.coordination_environments[0][0]['ce_symbol'][-1],
                    'ce_fraction':lse.coordination_environments[0][0]['ce_fraction'],
                    'csm':lse.coordination_environments[0][0]['csm'],
                    'permutation':lse.coordination_environments[0][0]['permutation']}
        return analysis
    
    def get_coordination(self, doc):
        """helper function that given an analyzed approx_neb doc, will extract and aggregate data into a
        dictionary with the keys as
        (1) input_csms  (2) input_coord_nums  (3) input_coord_envs
        (4) output_cms  (5) output_coord_nums (6) output_coord_envs
        and each value would be another dictionary with e.g. key = 0+1 and value = list of corresponding data
        """
        # initialize empty dictionaries
        input_csms, input_coord_nums, input_coord_envs, \
        output_csms, output_coord_nums, output_coord_envs = ({} for i in range(6))
        
        path_names = doc['images'].keys()  # get list of path names
        for path in path_names:
            # add end_points to images 
            start_point = doc['end_points'][int(path.split('+')[0])]
            final_point = doc['end_points'][int(path.split('+')[1])]
            path_docs = []
            path_docs.append(start_point)  # start point
            for d in doc['images'][path]:  # all the middle points
                path_docs.append(d)
            path_docs.append(final_point)  # final point
            
            # initialize empty lists to put aggregated data into
            input_csm, output_csm, input_cn, output_cn, input_ce, output_ce = ([] for i in range(6))
            
            # go through docs and sort into the corresponding lists
            for n, d in enumerate(path_docs):
                input_csm.append(get(d, ['input', 'csm']))
                output_csm.append(get(d, ['output', 'csm']))
                input_cn.append(get(d, ['input', 'coord_num']))
                output_cn.append(get(d, ['output', 'coord_num']))
                input_ce.append(get(d, ['input', 'coord_env']))
                output_ce.append(get(d, ['output', 'coord_env']))

            # add lists to dictionaries
            input_csms[path] = input_csm
            output_csms[path] = output_csm
            input_coord_nums[path] = input_cn
            output_coord_nums[path] = output_cn
            input_coord_envs[path] = input_ce
            output_coord_envs[path] = output_ce
            
        # add all the aggregated data into coordination dictionary for return
        coordination = {'input_csms':input_csms,
                       'output_csms':output_csms,
                       'input_coord_num':input_coord_nums,
                       'output_coord_num':output_coord_nums,
                       'intput_coord_envs':input_coord_envs,
                       'output_coord_envs':output_coord_envs}
        return coordination
    
    def unary_function(self, doc):
        """given a doc (dict) from the approx_neb collection, analyzes the coordination environment of 
        the sites along path (start site, path sites, end site) and returns the result in a new doc (dict)
        """
        wf_uuid = doc['wf_uuid']            # get wf_uuid for this doc
        obj_id = doc['_id']                 # get object id for this doc
        coord_env_analysis_output = {}      # initialize coord_env_analysis dictionary
        images_analyzed = {}                # initialize empty images dictionary
        print(f'\nAnalyzing doc (wf_uuid: {wf_uuid})')
    
        # get endpoints
        end_points = doc['end_points']  # list of stable sites
        
        # analyze endpoints
        if len(end_points) > 0:
            ion_index = doc['pathfinder'][list(doc['pathfinder'].keys())[0]]['relax_site_indexes'][0]  # find index of working ion
            elem = end_points[0]['input_structure']['sites'][ion_index]['label']   # find working ion
            end_points_analyzed = [None]*len(end_points)
            for n, site in enumerate(end_points):
                try:
                    input_struct = Structure.from_dict(site['input_structure'])         # get input structure
                    input_analyzed = self.analyze_struct_chemenv(input_struct, elem)    # analyze input structure
                except:
                    #print(f' Error: maybe missing input structure in endpoint {n}')
                    input_analyzed = None
                try:
                    output_struct = Structure.from_dict(site['output']['structure'])    # get output structure
                    output_analyzed = self.analyze_struct_chemenv(output_struct, elem)          # analyze output structure
                except:
                    #print(f' Error: maybe missing output structure in endpoint {n}')
                    output_analyzed = None
                end_points_analyzed[n] = {'input':input_analyzed, 'output':output_analyzed}
        else:
            end_points_analyzed = None

        # analyze images
        try:
            for path_name, path_sites in doc['images'].items():
                num_path_sites = len(path_sites)
                list_path_sites = [None] * num_path_sites
                for i in range(num_path_sites):
                    try:
                        input_struct = Structure.from_dict(path_sites[i]['input_structure'])      # get input structure
                        input_analyzed = self.analyze_struct_chemenv(input_struct, elem)                  # analyze input structure
                    except:
                        print(f' Error: maybe missing input structure in images {path_name} {i}')
                        input_analyzed = None
                    try:
                        output_struct = Structure.from_dict(path_sites[i]['output']['structure']) # get output structure
                        output_analyzed = self.analyze_struct_chemenv(output_struct, elem)                # analyze output structure
                    except:
                        print(f' Error: maybe missing output structure in images {path_name} {i}')
                        output_analyzed = None
                    list_path_sites[i] = {'input':input_analyzed, 'output':output_analyzed}
                images_analyzed[str(path_name)] = list_path_sites   # save analyzed path sites to path name
        except:
            print(f' Error analyzing images: potentially missing images to wf_uuid {wf_uuid}!')
            pass

        coord_env_analysis_output = {'images':images_analyzed, 'end_points':end_points_analyzed}
        
        #aggregating and sorting into coordination field
        coordination = self.get_coordination(coord_env_analysis_output)
        coord_env_analysis_output['coordination'] = coordination
        coord_env_analysis_output['wf_uuid'] = wf_uuid
        return coord_env_analysis_output
    
# query only specific types of doc
my_query = {'tags':{'$in':['20191122_batch']}}

chemenv_coord_env_builder = ChemEnvCoordEnvBuilder(approx_neb, sandbox, query=my_query, chunk_size=1)
chemenv_coord_env_builder.logger.disabled = True
chemenv_coord_env_builder.run()