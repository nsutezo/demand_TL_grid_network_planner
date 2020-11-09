import sys
import fiona
import os
import copy
import simone_agg_clustering as ac
import pandas as pd
import numpy as np
import string

specific_test_files ={('Kiriita','Nyandarua',3):2,
                      ('Ichagaki', 'Murang"A',3):7,
                      ('YimboWest','Siaya',3):0,
                      ('LusheyaLubinu','Kakamega',3):8,
                      ('AthiruRuujine', 'Meru',4):6,
                      ('AthiruRuujine', 'Meru',4):10,
                      ('Musanda', 'Kakamega',4):11,
                      ('Nyamaiya', 'Nyamira',4):12,
                      ('Githiga', 'Laikipia',3):7,
                      ('Masaba', 'Migori',2):0,
                      ('Kiguchwa', 'Meru',2):2
}


def count_structs_in_cell(src,bbox):
    count = 0
    for elem in src:
        if int(np.floor(elem['geometry']['coordinates'][0])) >=bbox[0] and int(np.floor(elem['geometry']['coordinates'][0] <bbox[2])):
            if int(np.floor(elem['geometry']['coordinates'][1])) >=bbox[1] and int(np.floor(elem['geometry']['coordinates'][1] <bbox[-1])):
                count+=1
    return count

def make_df(src,bbox):
    df =  []
    count = 0
    for elem in src:
        if int(np.floor(elem['geometry']['coordinates'][0])) >=bbox[0] and int(np.floor(elem['geometry']['coordinates'][0] <bbox[2])):
            if int(np.floor(elem['geometry']['coordinates'][1])) >=bbox[1] and int(np.floor(elem['geometry']['coordinates'][1] <bbox[-1])):
                df.append([elem['geometry']['coordinates'][0],elem['geometry']['coordinates'][1]])
                count+=1
    df = pd.DataFrame(df, columns=['x', 'y'])
    return df




def check_cell_size(src,scale):
    extent = src.bounds  #minx, miny, maxx, maxy
    y_step = (extent[-1] - extent[1])/float(scale)
    x_step = (extent[2] - extent[0])/float(scale)

    structs_counts_in_cell =[]
    for x in np.arange(extent[0], extent[2],x_step):
        for y in np.arange(extent[1], extent[-1], y_step):
            bbox = x,y,x+x_step,y+y_step
            bbox=[int(np.floor(b)) for b in bbox]
            structs_counts_in_cell.append(count_structs_in_cell(src,bbox))
    structs_counts_in_cell = np.array(structs_counts_in_cell)
    mean_count = np.mean(structs_counts_in_cell)
    #print(mean_count)
    if mean_count <= 1500:
        return True
    else:
        return False

def make_savetoken(cur_id, srcpath,cluster_idx):
    alphabet = list(string.ascii_lowercase)
    #subgrid = 'subgrid' +'_' +str(box[0]) +'_' +str(box[1]) +'_' +str(box[2]) +'_' +str(box[3])
    subgrid = 'subgrid_'+ str(cur_id)+'_'+ alphabet[cluster_idx]	
    if not os.path.exists(os.path.join(srcpath,subgrid)):
        os.mkdir(os.path.join(srcpath,subgrid))

    savetoken = os.path.join(srcpath,subgrid) +'/' + subgrid + '.shp'
    return savetoken

def make_small_shp(src,bbox,filename):
    output_schema = src.schema.copy()
    with fiona.open(filename, 'w', 'ESRI Shapefile', output_schema, crs=src.crs) as output:
        for elem in src:
            if int(np.floor(elem['geometry']['coordinates'][0])) >=bbox[0] and int(np.floor(elem['geometry']['coordinates'][0] <bbox[2])):
                if int(np.floor(elem['geometry']['coordinates'][1])) >=bbox[1] and int(np.floor(elem['geometry']['coordinates'][1] <bbox[-1])):
                    try:
                        output.write(elem)
                    except Exception as err:
                        print('Error writing following',elem)

def split_iterative(start_bbox,src):
    stack, valid = [start_bbox],[]
    while stack:
        bbox = stack.pop()
        print(bbox)
        cur_df = make_df(src,bbox)
	try:
            max_nodes = ac.get_max_nodes_in_cluster(cur_df)
            print('Max Nodes: ', max_nodes)
            if max_nodes < 250:
                valid.append(bbox)
            else:
                y_step = (bbox[-1] - bbox[1])/2.0
                x_step = (bbox[2] - bbox[0])/2.0
                for x in np.arange(bbox[0], bbox[2],x_step):
                    for y in np.arange(bbox[1], bbox[-1], y_step):
                        sub_bbox = x,y,x+x_step,y+y_step
                        sub_bbox=[int(np.floor(b)) for b in sub_bbox]
                        stack.append(sub_bbox)
                        print('Stack size: ',len(stack))
        except Exception as err:
	    print(err)
    return valid

if __name__ =='__main__':
    path ='../../selected_wards_gypsum/'
    save_path = '../../iterative_split/'
    wards_in_gypsum = os.listdir(path)
    wards_in_gypsum = [(r.split('_')[0],r.split('_')[1]) for r in wards_in_gypsum]
    #for w_c_s,subgrid in specific_test_files.iteritems():
    #w_c_s=sys.argv[1]
    #print('Starting with:', w_c_s)
    #import ipdb; ipdb.set_trace()
    w,c,s = sys.argv[1], sys.argv[2], int(sys.argv[3])
    w_c_s = (w,c,s)
    subgrid = int(sys.argv[4])
    print('Ward ,County, Scale',w,c,s)
    cur_path= path+ w+'_'+c +'/scale_' + str(s) + '/subgrid_' + str(subgrid) + '/'
    shp_path = cur_path + 'subgrid_' + str(subgrid) +'.shp'
    save_path = save_path+ w+'_'+c
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    save_path += '/scale_' + str(s)
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    save_path += '/subgrid_' + str(subgrid)
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    save_path+='/'
    src = fiona.open(shp_path)
    extent = src.bounds
    y_step = (extent[-1] - extent[1])/2.0
    x_step = (extent[2] - extent[0])/2.0
    idx =0
    for x in np.arange(extent[0], extent[2],x_step):
        for y in np.arange(extent[1], extent[-1], y_step):
            sub_bbox = x,y,x+x_step,y+y_step
            sub_bbox=[int(np.floor(b)) for b in sub_bbox]
            print('Current subbox ',sub_bbox)
            print('Running iterative split')
            valid_futher_splits = split_iterative(sub_bbox,src)
            print('Sample acceptable extent',valid_futher_splits[0])
            print('Starting the make shapefile process')
            for further_split in valid_futher_splits:
                cur_df = make_df(src,further_split)
                savetoken = make_savetoken(subgrid,save_path,idx)
                make_small_shp(src,further_split,savetoken)
                print('Save shapefile for',w_c_s, further_split)
		idx+=1

print('Done')

