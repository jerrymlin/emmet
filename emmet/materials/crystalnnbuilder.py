import json
import os
import urllib as url
from pprint import pprint
from pymatgen import Structure
from pymatgen.analysis.local_env import CrystalNN
from pydash import get
from pymongo import MongoClient

from tqdm import tqdm
from maggma.core import Builder
from maggma.core import Store
from maggma.builders.map_builder import MapBuilder
from maggma.stores import MongograntStore

# Connecting to approx_neb database
approx_neb = MongograntStore("ro:mongodb03.nersc.gov/fw_acr_mv","approx_neb")
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

class CoordEnvBuilder(MapBuilder):
    """
    given a doc (dict) from the approx_neb collection, analyzes the coordination environment of 
    the sites along path (start site, path sites, end site) and returns the result in a new doc (dict)
    
    Args:
        source_collection[MongoStore]: collection specified by MongoStore() from maggma.stores, input collection
        target_collection[MongoStore]: collection specified by MongoStore() from maggma.stores, output collection
    """
    def get_index_and_coords(self, struct, elem):
        """helper function that gets the coordinates and index of a given element in a given structure"""
        for n, site in enumerate(struct.sites):
            if str(site.specie) == elem:
                index = n
                coords = struct.sites[n].coords
        # safety check for if working ion was found in structure sites
        try:
            index
        except NameError:
            print('Could not find working ion in structure sites!')
        return index, coords
    
    def analyze_struct(self, s, elem):
        """helper function that given a structure (s) and a working ion (elem), analyzes the coordination 
        environment and returns the output as a dict"""
        cnn = CrystalNN()
        index, coords = self.get_index_and_coords(s, elem)    # find index (and coordinates) of element
        cn = cnn.get_cn(s, index)                        # get coordination number w/ index
        op = cnn.get_local_order_parameters(s, index)    # get order parameters w/ index
        max_op = max(op.values())                        # determine structure motif
        for struct_type, order_param in op.items():
            if order_param == max_op:
                found_struct_type = struct_type

        analysis = {'structure':s.as_dict(),               # save results into output dictionary
                    'working_ion':elem,
                    'order_parameters':op,
                    'coord_num':cn,
                    'structure_motif':found_struct_type}
        return analysis
    
    def get_coordination(self, doc):
        """helper function that given an analyzed approx_neb doc, will extract and aggregate data into a
        dictionary with the keys as
        (1) input_order_parameters  (2) input_coord_num  (3) input_structure_motif
        (4) output_order_parameters (5) output_coord_num (6) output_structure_motif
        and each value would be another dictionary with e.g. key = 0+1 and value = list of corresponding data
        """
        # initialize empty dictionaries
        input_order_parameters, input_coord_num, input_structure_motif, \
        output_order_parameters, output_coord_num, output_structure_motif = ({} for i in range(6))
        
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
            input_op, output_op, input_cn, output_cn, input_sm, output_sm = ([] for i in range(6))
            
            # go through docs and sort into the corresponding lists
            for n, d in enumerate(path_docs):
                input_op.append(get(d, ['input', 'order_parameters']))
                output_op.append(get(d, ['output', 'order_parameters']))
                input_cn.append(get(d, ['input', 'coord_num']))
                output_cn.append(get(d, ['output', 'coord_num']))
                input_sm.append(get(d, ['input', 'structure_motif']))
                output_sm.append(get(d, ['output', 'structure_motif']))

            # add lists to dictionaries
            input_order_parameters[path] = input_op
            output_order_parameters[path] = output_op
            input_coord_num[path] = input_cn
            output_coord_num[path] = output_cn
            input_structure_motif[path] = input_sm
            output_structure_motif[path] = output_sm
            
        # add all the aggregated data into coordination dictionary for return
        coordination = {'input_order_parameters':input_order_parameters,
                       'output_order_parameters':output_order_parameters,
                       'input_coord_num':input_coord_num,
                       'output_coord_num':output_coord_num,
                       'intput_structure_motif':input_structure_motif,
                       'output_structure_motif':output_structure_motif}
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
        try:
            end_points = doc['stable_sites']  # list of stable sites
        except:
            end_points = doc['end_points']    # list of end_points
        
        # analyze endpoints
        if len(end_points) > 0:
            
            ion_index = doc['pathfinder'][list(doc['pathfinder'].keys())[0]]['relax_site_indexes'][0]  # find index of working ion
            elem = end_points[0]['input_structure']['sites'][ion_index]['label']   # find working ion
            end_points_analyzed = [None]*len(end_points)
            for n, site in enumerate(end_points):
                try:
                    input_struct = Structure.from_dict(site['input_structure'])         # get input structure
                    input_analyzed = self.analyze_struct(input_struct, elem)            # analyze input structure
                except:
                    #print(f' Error: maybe missing input structure in endpoint {n}')
                    input_analyzed = None
                try:
                    output_struct = Structure.from_dict(site['output']['structure'])    # get output structure
                    output_analyzed = self.analyze_struct(output_struct, elem)          # analyze output structure
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
                        input_analyzed = self.analyze_struct(input_struct, elem)                  # analyze input structure
                    except:
                        #print(f' Error: maybe missing input structure in images {path_name} {i}')
                        input_analyzed = None
                    try:
                        output_struct = Structure.from_dict(path_sites[i]['output']['structure']) # get output structure
                        output_analyzed = self.analyze_struct(output_struct, elem)                # analyze output structure
                    except:
                        #print(f' Error: maybe missing output structure in images {path_name} {i}')
                        output_analyzed = None
                    list_path_sites[i] = {'input':input_analyzed, 'output':output_analyzed}
                images_analyzed[str(path_name)] = list_path_sites   # save analyzed path sites to path name
        except:
            #print(f' Error analyzing images: potentially missing images to wf_uuid {wf_uuid}!')
            pass

        coord_env_analysis_output = {'images':images_analyzed, 'end_points':end_points_analyzed}
        
        #aggregating and sorting into coordination field
        coordination = self.get_coordination(coord_env_analysis_output)
        coord_env_analysis_output['coordination'] = coordination
        coord_env_analysis_output['wf_uuid'] = wf_uuid
        return coord_env_analysis_output
    
# query only specific types of doc
my_query = {'batt_id':{'$regex':'Mg'},
            'tags':{'$in':['20191122_batch']}, 
            'last_updated':{'$exists':1}}

coord_env_builder = CoordEnvBuilder(source_collection, target_collection, query=my_query)

# run builder
coord_env_builder.run()