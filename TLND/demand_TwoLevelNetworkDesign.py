#!/usr/bin/env python
"""
A heuristic approach for two-level network design - rural electrification
Ayse Selin Kocaman
ask2170@columbia.edu
"""

import os, time, copy
import CMST_dfs_OLD
import gc
import collections
from heapq import heappush, heappop
from osgeo import ogr
import network
import fileRW
import numpy as np
import shutil
# import fiona
import pandas as pd


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


class Error(Exception):
    def __init__(self, msg):
        self.msg = msg


def mergeCluster(ClusterByNode, NodesByClusterID, Centers, segment):
    center1, center2 = segment.getNodes()

    centerX = (center1.getWeight() * center1.getCenterX()
               + center2.getWeight() * center2.getCenterX()) / (center2.getWeight() + center1.getWeight())
    centerY = (center1.getWeight() * center1.getCenterY()
               + center2.getWeight() * center2.getCenterY()) / (center2.getWeight() + center1.getWeight())

    weight = center2.getWeight() + center1.getWeight()
    baseClusterID = min(ClusterByNode[center1], ClusterByNode[center2])
    mergingClusterID = max(ClusterByNode[center1], ClusterByNode[center2])

    NodesByClusterID[baseClusterID].extend(NodesByClusterID.pop(mergingClusterID))

    Centers[baseClusterID].setXY(centerX, centerY)
    Centers[baseClusterID].setWeight(weight)

    del Centers[mergingClusterID]

    for node in NodesByClusterID[baseClusterID]:
        ClusterByNode[node] = baseClusterID


def generateDictsFromShp(shapeFile, outputPath):
    'Reads nodes and node weights from a point shapefile.'
    rootDir, fc = os.path.split(shapeFile)
    file, ext = os.path.splitext(fc)

    if not os.path.exists(outputPath):
        try:
            os.mkdir(outputPath)
        except:
            print("ERROR: could not create new directory", outputPath)
    ds = ogr.Open(shapeFile)
    ptLayer = ds.GetLayer(0)

    nodesByClusterID = collections.defaultdict(list)
    clusterByNode = {}
    nodes = {}
    centers = {}
    LVCostDict = {}
    # X = np.array([[]]).reshape(0,4)
    feat = ptLayer.GetNextFeature()
    #    nodes_weights_output = []
    nodes_demands_output = []
    # import ipdb; ipdb.set_trace()
    np.random.seed(7)
    indices = np.random.permutation(ptLayer.GetFeatureCount())
    high = indices[:int(0.3 * ptLayer.GetFeatureCount())]

    while feat is not None:
        FID = feat.GetFID()
        nodeWeight = 1  # np.random.randint(low=1,high=101)# testing the effect of weights on the design 1
        if FID in high:
            nodeDemand = 100  # at 30 kWh/ month
        else:
            nodeDemand = 30  # np.random.randint(low=60,high=150) # at 60 kWh/month
        geomRef = feat.GetGeometryRef()
        x = geomRef.GetX()
        y = geomRef.GetY()
        FID = feat.GetFID()
        nodes[FID] = network.Node(FID, x, y, nodeWeight, nodeDemand)  # Households
        centers[FID] = network.Node(FID, x, y, nodeWeight, nodeDemand)  # Transformers (center of mass of the cluster)
        # X =np.concatenate((X,np.array([[FID,x,y,nodeWeight]])), axis = 0)
        clusterByNode[nodes[FID]] = FID
        nodesByClusterID[FID].append(nodes[FID])
        LVCostDict[FID] = 0
        #        nodes_weights_output.append([x,y,nodeWeight])
        nodes_demands_output.append([x, y, nodeDemand])
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    nodes_demands_output = pd.DataFrame(nodes_demands_output, columns=['x', 'y', 'demands'])
    return nodesByClusterID, clusterByNode, nodes, centers, LVCostDict, nodes_demands_output


def generateSegments(centers, searchRadius):
    segments = []
    nodeCopy = centers.copy()

    segID = 0
    for startNode in centers.values():
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            dist = ((startNode.getX() - endNode.getX()) ** 2 +
                    (startNode.getY() - endNode.getY()) ** 2) ** (.5)
            if dist < searchRadius:
                segments.append(
                    network.Seg(segID, startNode, endNode, dist, 1))  # 1 at end because added a param for Demand
                segID += 1
    return segments


def generateSegmentsDemand(centers, searchRadius, all_ready_checked_list, roi_years, cost_per_kwh, tcost, LV,
                           fraction_recovered):
    segments = []
    nodeCopy = centers.copy()
    # multiplier = 1
    segID = 0

    for startNode in centers.values():
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            # if (startNode.getID() == 605) and (endNode.getID() == 606):
            # import ipdb; ipdb.set_trace()
            #            if (startNode.getID() !=node1) & (endNode.getID()!=node2):

            if [startNode.getID(), endNode.getID()] not in all_ready_checked_list:
                dist = ((startNode.getX() - endNode.getX()) ** 2 +
                        (startNode.getY() - endNode.getY()) ** 2) ** (.5)
                total_demand = startNode.getDemand() + endNode.getDemand()  # calculates the demand from both nodes
                revenue = total_demand * roi_years * 12 * cost_per_kwh
                # tx_cost = tcost * ((total_demand * 1000 * 0.6) / (4.0 * 30.0)) * fraction_recovered
                delta = revenue - tcost
                if delta > 0:
                    searchRadius = delta / float(LV * fraction_recovered)
                    if dist < searchRadius:  # selects feasible segments
                        segments.append(network.Seg(segID, startNode, endNode, dist, total_demand))
                        segID += 1
            else:
                pass
    #                print [startNode.getID(),endNode.getID()]
    # check segments length
    #    while (len(segments) == 0) & (multiplier <=5) :
    #        segments = []
    #        nodeCopy = centers.copy()
    #        segID = 0
    #        for startNode in centers.values():
    #            del nodeCopy[startNode.getID()]
    #            for endNode in nodeCopy.values():
    #                dist = ((startNode.getX() - endNode.getX()) ** 2 +
    #                        (startNode.getY() - endNode.getY()) ** 2) ** (.5)
    #                total_demand =  startNode.getDemand() + endNode.getDemand() # calculates the demand from both nodes
    #                if dist < searchRadius*multiplier: #selects feasible segments
    #                    segments.append(network.Seg(segID, startNode, endNode,dist,total_demand))
    #                    segID += 1
    #        multiplier+=1

    return segments


def generateSegmentsDemand_debug(centers, searchRadius, all_ready_checked_list):
    segments = []
    nodeCopy = centers.copy()
    multiplier = 1
    segID = 0

    for startNode in centers.values():
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            # if (startNode.getID() == 605) and (endNode.getID() == 606):
            # import ipdb; ipdb.set_trace()
            #            if (startNode.getID() !=node1) & (endNode.getID()!=node2):

            if [startNode.getID(), endNode.getID()] not in all_ready_checked_list:
                dist = ((startNode.getX() - endNode.getX()) ** 2 +
                        (startNode.getY() - endNode.getY()) ** 2) ** (.5)
                total_demand = startNode.getDemand() + endNode.getDemand()  # calculates the demand from both nodes
                if dist < searchRadius:  # selects feasible segments
                    segments.append(network.Seg(segID, startNode, endNode, dist, total_demand))
                    segID += 1
            else:
                pass
    return segments


def generateSegmentsDemand2(centers, searchRadius, all_ready_checked_list_mv):
    segments = []
    nodeCopy = centers.copy()
    segID = 0
    # all_ready_checked_list_mv = set(all_ready_checked_list_mv)
    for startNode in centers.values():
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            # if (startNode.getID() == 680)  and (endNode.getID() ==683):
            # import ipdb; ipdb.set_trace()
            if [startNode.getID(), endNode.getID()] not in all_ready_checked_list_mv:
                dist = ((startNode.getX() - endNode.getX()) ** 2 +
                        (startNode.getY() - endNode.getY()) ** 2) ** (.5)
                total_demand = startNode.getDemand() + endNode.getDemand()  # calculates the demand from both nodes
                if (dist < searchRadius):  # selects feasible segments
                    segments.append(network.Seg(segID, startNode, endNode, dist, total_demand))
                    segID += 1
            else:
                pass
    return segments


def maxInClusterDist(centerNode, nodesByClusterID):  # Returns maxDist within the cluster
    maxdist = 0
    for node in nodesByClusterID[centerNode.getID()]:  # uses the fact that centerID and ClusterID are same
        dist = ((centerNode.getX() - node.getX()) ** 2 +
                (centerNode.getY() - node.getY()) ** 2) ** (.5)
        if dist >= maxdist:
            maxdist = dist
    return maxdist


def maxTempInClusterDist(segment, ClusterByNode, nodesByClusterID, distFromT):
    flag = True
    tempCenter1, tempCenter2 = segment.getNodes()

    tempCenterX = (tempCenter1.getWeight() * tempCenter1.getX()
                   + tempCenter2.getWeight() * tempCenter2.getX()) / (tempCenter2.getWeight() + tempCenter1.getWeight())
    tempCenterY = (tempCenter1.getWeight() * tempCenter1.getY()
                   + tempCenter2.getWeight() * tempCenter2.getY()) / (tempCenter2.getWeight() + tempCenter1.getWeight())
    # given the cluster ID of the node, find all nodes in the cluster.
    # for each node in cluster calculate the distance between centroid and the node
    # Find the largest distance
    for node in nodesByClusterID[ClusterByNode[segment.getNode1()]]:
        dist = ((tempCenterX - node.getX()) ** 2 + (tempCenterY - node.getY()) ** 2) ** (.5)
        #        maxDemand += node.getWeight()
        if dist > distFromT:
            flag = False
            break
    if flag:
        for node in nodesByClusterID[ClusterByNode[segment.getNode2()]]:
            dist = ((tempCenterX - node.getX()) ** 2 + (tempCenterY - node.getY()) ** 2) ** (.5)
            if dist > distFromT:
                flag = False
                break

    return flag, tempCenterX, tempCenterY  # maxDist, maxDemand


def loggers(log_filename, initial=False, data=None):
    if initial:
        with open(log_filename, 'w') as src:
            src.write("")
    else:
        with open(log_filename, 'a') as src:
            src.write(str(data) + "\n")


def totalInClusterCost(nodesByClusterID, centers):
    totalCost = 0
    for centerID in centers.keys():
        for node in nodesByClusterID[centerID]:
            totalCost += ((node.getX() - centers[centerID].getX()) ** 2 +
                          (node.getY() - centers[centerID].getY()) ** 2) ** (.5)
    return totalCost


def kruskalsAlg(segments, nodes):
    'Kruskal\'s algorithm for finding a minimum spanning tree'
    segments.sort(key=lambda obj: obj.getWeight())
    tree = network.Network()
    numNodes = len(nodes)

    for segment in segments:
        node1 = segment.getNode1()
        node2 = segment.getNode2()
        node1InNet = tree.inNet(node1)
        node2InNet = tree.inNet(node2)

        if (not node1InNet and not node2InNet) or (node1InNet != node2InNet):
            tree.addSeg(segment)
        else:
            if node1InNet and node2InNet and \
                    (tree.getNetID(node1) != tree.getNetID(node2)):
                tree.addSeg(segment)
        if tree.numNodes() > numNodes:
            break
    return tree, segments


def primsAlg(numNodes, firstNodeID, nodeDict):
    'Prim\'s Algorithm for finding a minimum spanning tree'

    tree = network.Network()
    segHeap = []

    # Find the shortest segment emanating from the node with the firstNodeID
    try:
        segs = nodeDict[firstNodeID]
    except KeyError:
        return tree

    leastWeight = None
    for seg in segs:
        if (leastWeight == None):
            leastWeight = seg.getWeight()
            firstSeg = seg
        elif (seg.getWeight() < leastWeight):
            leastWeight = seg.getWeight()
            firstSeg = seg
    tree.addSeg(firstSeg)

    # Starter to algorithm
    # Add the segs emanating from the first two endpoints to the heap
    for endNode in [firstSeg.getNode1(), firstSeg.getNode2()]:
        addToHeap(segHeap, nodeDict[endNode.getID()])

    # Pick best from heap and repeat
    while tree.numNodes() < numNodes:
        try:
            # Get best segment from heap
            seg = heappop(segHeap)
        except:
            # Tree is finished (not all nodes contained).
            break
        node1, node2 = seg.getNodes()
        node1InNet = tree.inNet(node1)
        node2InNet = tree.inNet(node2)
        # Add the seg if it's terminal node isn't already in the cluster.
        if (not node1InNet) or (not node2InNet):
            if not node1InNet:
                endNode = node1
            else:
                endNode = node2
            tree.addSeg(seg)
            # Add all emanating segs to the heap:
            # nodeDict returns all segments coming out from the endNode
            # endNode is the node that is outside of the tree
            addToHeap(segHeap, nodeDict[endNode.getID()])
            # And we are sure that everything in the heap is adjacent to the tree because
            # we only add the adjacent segments in the first place using nodeDict
    return tree


def addToHeap(heap, newSegs):
    'Adds new segments to the segHeap.'
    for seg in newSegs:
        heappush(heap, seg)
    return heap


def buildAssocDict(segments):
    'Builds dictionary with nodeID key where values are segs from/to that node'
    segList = {}
    for seg in segments:
        node1, node2 = seg.getNodes()
        for nodeID in [node1.getID(), node2.getID()]:
            if nodeID in segList.keys():
                # if segList.has_key(nodeID):
                segList[nodeID].append(seg)
            else:
                segList[nodeID] = [seg]
    return segList


def convertDictToArray(centers):
    # import ipdb
    # ipdb.set_trace()
    # array = np.array([[]]).reshape(0,2)
    array = np.zeros((len(centers), 2))
    for c in range(len(centers)):
        node = centers[c]
        # array = np.concatenate((array,np.array([[node.getX(),node.getY()]])), axis = 0)
        # array[c,0],array[c,1] =node.getX(),node.getY()
        array[c, :] = [node.getX(), node.getY()]
    return array


def get_txcost(centers, nodesByClusterID, tcost):
    number_transformers = 0
    number_connected = 0
    connected_centers = {}
    tx_cost = 0
    connected_demand = 0
    for i, k in enumerate(centers):
        # check if the center has more than one node
        if len(nodesByClusterID[k]) > 1:
            number_connected += len(nodesByClusterID[k])
            number_transformers += 1
            tx_cost += tcost
            connected_centers[k] = centers[k]
            connected_demand += centers[k].getDemand()
    transformer_cost = tcost * ((connected_demand * 1000 * 0.6) / (4.0 * 30.0))  # assuming 4 hours of peak
    return number_transformers, number_connected, connected_centers, transformer_cost, connected_demand


def get_connected_nodes(txs, nodesByClusterID):
    connected_nodes_by_id = {}
    for i, k in enumerate(txs):
        if len(nodesByClusterID[k]) > 1:
            connected_nodes_by_id[k] = nodesByClusterID[k]
    return connected_nodes_by_id


def get_all_txcost(centers, nodesByClusterID, tcost, number_in_file):
    number_transformers = 0
    number_connected = 0
    connected_centers = {}
    tx_cost = 0
    connected_demand = 0
    for i, k in enumerate(centers):
        # check if the center has more than one node
        #        if len(nodesByClusterID[k])>1:
        number_connected += len(nodesByClusterID[k])
        number_transformers += 1
        tx_cost += tcost / float(centers[k].getDemand() * len(nodesByClusterID[k]))  # cost per KW per HH
        connected_centers[k] = centers[k]
        connected_demand += centers[k].getDemand()
    #    import ipdb;ipdb.set_trace()
    return tx_cost * connected_demand * number_in_file


def get_distance(centroidX, centroidY, nodelist):
    dist = 0
    for node in nodelist:
        tmp_dist = ((centroidX - node.getX()) ** 2 + (centroidY - node.getY()) ** 2) ** (.5)
        dist += tmp_dist
    return dist


def cost_converged(cost_list, threshold=0.01):
    for j in range(1, len(cost_list)):
        diff = abs(cost_list[j] - cost_list[j - 1]) / cost_list[j]
        if diff <= threshold:
            convergence = True
        else:
            convergence = False
            break
    return convergence


def get_worse_case_prims(tx_locations, MV, sr):
    try:
        # import ipdb; ipdb.set_trace()
        segments_ST = generateSegmentsDemand2(tx_locations, sr*1000000000, [])  #
        nodeDict_ST = buildAssocDict(segments_ST)
        minTree = primsAlg(len(tx_locations), [*nodeDict_ST.keys()][0], nodeDict_ST)
        minTotalCost_ST = minTree.getTotalEdgeWeight() * MV

    except:
        minTotalCost_ST = 100000000000
        minTree = None

    return minTotalCost_ST, minTree


def get_smallest_possible_MV_cost(tx_locations, MV, LV, minNodesByClusterID, minClusterByNode, distFromT, sr):
    # import ipdb; ipdb.set_trace()

    try:
        segments_ST = generateSegmentsDemand2(tx_locations, sr, [])  #
        nodeDict_ST = buildAssocDict(segments_ST)
        minTree_ST = primsAlg(len(tx_locations), [*nodeDict_ST.keys()][0], nodeDict_ST)
    except:
        minTotalCost_ST = 100000000000
        return minTotalCost_ST, None

    centers_ST = copy.deepcopy(tx_locations)
    minTotalCost_ST = minTree_ST.getTotalEdgeWeight() * MV

    nodesByClusterID_ST = copy.deepcopy(minNodesByClusterID)
    clusterByNode_ST = copy.deepcopy(minClusterByNode)

    try:
        minSeg_ST = min(segments_ST, key=lambda obj: obj.getWeight())
        flag, tempCenterX, tempCenterY = maxTempInClusterDist(minSeg_ST, clusterByNode_ST, nodesByClusterID_ST,
                                                              distFromT)
        if not flag:
            segments_ST.sort(key=lambda obj: obj.getWeight())

            for seg in segments_ST:
                if seg.getWeight() > distFromT * 2:
                    break  # break from for loop
                else:
                    flag, tempCenterX, tempCenterY = maxTempInClusterDist(seg, clusterByNode_ST, nodesByClusterID_ST,
                                                                          distFromT)
                    if flag:
                        minSeg_ST = seg
                        break  # break from for loop
    except:
        return minTotalCost_ST

    if minSeg_ST.getWeight() <= distFromT * 2:
        flag = True  # can be anything less than 500
    else:
        flag = False  # can be anything greater than 500
        print("NO CLUSTER POSSIBLE")

    all_ready_checked_mv = []

    tempCenter1, tempCenter2 = minSeg_ST.getNodes()

    tempCenterX = (tempCenter1.getWeight() * tempCenter1.getX()
                   + tempCenter2.getWeight() * tempCenter2.getX()) / (tempCenter2.getWeight() + tempCenter1.getWeight())
    tempCenterY = (tempCenter1.getWeight() * tempCenter1.getY()
                   + tempCenter2.getWeight() * tempCenter2.getY()) / (tempCenter2.getWeight() + tempCenter1.getWeight())

    counter = len(tx_locations)

    delta = 100
    threshold = 0.1
    cost_flag = minTotalCost_ST
    while flag and delta > threshold:  # and (minTotalCost_ST > wiggle_room):

        tmp_nodesByClusterID_ST = copy.deepcopy(nodesByClusterID_ST)
        tmp_centers_ST = copy.deepcopy(centers_ST)
        tmp_clusterByNode_ST = copy.deepcopy(clusterByNode_ST)
        # tmp_LVCostDict = copy.deepcopy(LVCostDict)

        center1, center2 = minSeg_ST.getNodes()

        weight = center2.getWeight() + center1.getWeight()
        baseClusterID_ST = min(clusterByNode_ST[center1], clusterByNode_ST[center2])

        mergingClusterID_ST = max(clusterByNode_ST[center1], clusterByNode_ST[center2])
        nodesByClusterID_ST[baseClusterID_ST].extend(nodesByClusterID_ST.pop(mergingClusterID_ST))

        centers_ST[baseClusterID_ST].setXY(tempCenterX, tempCenterY)
        centers_ST[baseClusterID_ST].setWeight(weight)

        all_ready_checked_mv.append([mergingClusterID_ST, baseClusterID_ST])
        del centers_ST[mergingClusterID_ST]

        for node in nodesByClusterID_ST[baseClusterID_ST]:
            clusterByNode_ST[node] = baseClusterID_ST

        # generate segments for new graph
        segments_ST = generateSegmentsDemand2(centers_ST, sr, all_ready_checked_mv)
        # if len(segments_ST)>0:
        nodeDict_ST = buildAssocDict(segments_ST)
        newTree_ST = primsAlg(len(centers_ST), [*nodeDict_ST.keys()][0],
                              nodeDict_ST)  # 0 is the starting node of prims.
        TotalMVCost_ST = newTree_ST.getTotalEdgeWeight() * MV
        print('MV cost:', TotalMVCost_ST)

        gc.collect()

        additional_LV = get_distance(tempCenterX, tempCenterY, [center1, center2])
        newTotalCost_ST = TotalMVCost_ST + additional_LV * LV
        print(counter)
        counter -= 1

        delta = abs(cost_flag - newTotalCost_ST)  # / minTotalCost_ST
        cost_flag = newTotalCost_ST

        print('All ready checked segments:', len(all_ready_checked_mv))
        print('New cost:', newTotalCost_ST)
        if (newTotalCost_ST <= minTotalCost_ST):
            minTotalCost_ST = newTotalCost_ST
            minTree_ST = copy.deepcopy(newTree_ST)
        else:
            centers_ST = copy.deepcopy(tmp_centers_ST)
            nodesByClusterID_ST = copy.deepcopy(tmp_nodesByClusterID_ST)
            clusterByNode_ST = copy.deepcopy(tmp_clusterByNode_ST)
            segments_ST = generateSegmentsDemand_debug(centers_ST, sr, all_ready_checked_mv)
        print('Number of Segments:', len(segments_ST))
        # Calculate maxDist below for next graph and continue if it is less than 500

        try:  # to check if there is a segment on the graph or there is only one cluster
            minSeg_ST = min(segments_ST, key=lambda obj: obj.getWeight())
            flag, tempCenterX, tempCenterY = maxTempInClusterDist(minSeg_ST, clusterByNode_ST, nodesByClusterID_ST,
                                                                  distFromT)
            if not flag:
                segments_ST.sort(key=lambda obj: obj.getWeight())

                for seg in segments_ST:
                    if seg.getWeight() > distFromT * 2:
                        break
                    else:
                        flag, tempCenterX, tempCenterY = maxTempInClusterDist(seg, clusterByNode_ST,
                                                                              nodesByClusterID_ST,
                                                                              distFromT)
                        if flag:
                            minSeg_ST = seg
                            break
        except Exception as err:
            # print(err)
            break

    return minTotalCost_ST, minTree_ST, centers_ST, nodesByClusterID_ST


def run(centers, nodesByClusterID, clusterByNode, LVCostDict, sr, MV, LV, TCost, distFromT, investment, cost_per_kwh,
        roi_years, fraction_recovered, maxLVLenghtInCluster,
        outputDir, logfilename):
    print("JOINT LV , MV optimization")
    minCenters = copy.deepcopy(centers)
    # fraction_lv_recovered = 0.5
    # fraction_mv_recovered = 1 - fraction_lv_recovered
    # fraction_lv_investment = 0.5
    sumLVCostAtEachStep = {}
    segments = generateSegmentsDemand(minCenters, sr, [], roi_years, cost_per_kwh, TCost, LV,
                                      fraction_recovered)  # find closest points within 1000m radius
    minSeg = min(segments, key=lambda obj: obj.getWeight())  # find the 2 closest nodes
    # To write total cost to a text file
    statFile = outputDir + os.sep + "TotalCost_FirstStage.txt"
    outFile = open(statFile, "w")

    if minSeg.getWeight() <= distFromT * 2:  # if dist between 2 closest nodes <1000
        # maxDist = 0  # can be anything less than 500
        flag = True
    else:
        # maxDist = maxLVLenghtInCluster + 10  # can be anything greater than 500
        flag = False

    # revenue = minSeg.getDemand() * roi_years * 12 * cost_per_kwh

    #    print "NO CLUSTER POSSIBLE"
    # find the centroid of 2 closest points
    tempCenter1, tempCenter2 = minSeg.getNodes()
    tempCenterX = (tempCenter1.getWeight() * tempCenter1.getX()
                   + tempCenter2.getWeight() * tempCenter2.getX()) / (tempCenter2.getWeight() + tempCenter1.getWeight())
    tempCenterY = (tempCenter1.getWeight() * tempCenter1.getY()
                   + tempCenter2.getWeight() * tempCenter2.getY()) / (tempCenter2.getWeight() + tempCenter1.getWeight())
    i = len(centers)
    initial = True
    loggers(logfilename, initial)
    initial = False
    # num_tx = 0
    total_cost = 0

    all_ready_checked = []
    while flag:  # & (total_cost < investment):  # (maxDist <= distFromT):
        #    while (maxDist <= distFromT) & (newTotalCost <= investment):
        i -= 1  # update the value of i

        if i % 20 == 0:
            print(i)
        # if i == 698:
        # import ipdb; ipdb.set_trace()

        cur_token = 'stage1 ' + str(i)
        loggers(logfilename, initial, cur_token)
        tmp_nodesByClusterID = copy.deepcopy(nodesByClusterID)
        tmp_centers = copy.deepcopy(centers)
        tmp_clusterByNode = copy.deepcopy(clusterByNode)
        tmp_LVCostDict = copy.deepcopy(LVCostDict)

        center1, center2 = minSeg.getNodes()  # get 2 closest nodes
        weight = center2.getWeight() + center1.getWeight()  # get the weights
        demand = center2.getDemand() + center1.getDemand()
        baseClusterID = min(clusterByNode[center1], clusterByNode[center2])  # picks smallest cluster id
        mergingClusterID = max(clusterByNode[center1], clusterByNode[center2])  # picks the largest cluster id
        nodesByClusterID[baseClusterID].extend(
            nodesByClusterID.pop(mergingClusterID))  # merges the larger id cluster to the smaller cluster id.
        centers[baseClusterID].setXY(tempCenterX, tempCenterY)  # sets centroid as the new xy for the smaller cluster id
        centers[baseClusterID].setWeight(weight)  # sets the new weight
        centers[baseClusterID].setDemand(demand)
        #        print mergingClusterID
        del centers[mergingClusterID]  # removes merged node's id from centers

        for node in nodesByClusterID[baseClusterID]:  # assigns the same cluster id to every node in clusterByNode
            clusterByNode[node] = baseClusterID

        #        TotalTransformerCost = len(centers) * TCost
        number_transformers, number_connected, tx_locations, TotalTransformerCost, connected_demand = get_txcost(
            centers, nodesByClusterID, TCost)  # updates the tx cost
        #        TotalTransformerCost = number_transformers * TCost

        del LVCostDict[mergingClusterID]  # removes the LV cost of the merged cluster
        gc.collect()
        # given the new centroids, it calculates the LV params
        segmentsCMST, LVCostDict[baseClusterID] = CMST_dfs_OLD.CMST(nodesByClusterID[baseClusterID],
                                                                    maxLVLenghtInCluster,
                                                                    centers[baseClusterID])

        sumLVCostAtEachStep[len(centers)] = sum(LVCostDict.values()) * LV
        revenue = connected_demand * roi_years * 12 * cost_per_kwh  # * fraction_lv_recovered# dollar amount/per day

        # mv_wiggle_room = (revenue - (TotalTransformerCost + (sum(LVCostDict.values())) * LV ) * fraction_recovered)/fraction_recovered
        # mv_wiggle_room = (revenue - (TotalTransformerCost + (sum(LVCostDict.values())) * LV))
        # Here we get the least cost MV (given that some tx can be regrouped).
        # If we can recover cost under the least cost MV then we accept merge
        # Knowing that at the end we can run the get_smallest_possible_MV to obtain MV
        # if len(tx_locations) > 1:
        #     lowest_MVCost, mv_tree = get_worse_case_prims(tx_locations, MV, sr+10)
        # lowest_MVCost, mv_tree = get_smallest_possible_MV_cost(tx_locations, MV,LV,
        #                                                nodesByClusterID, clusterByNode,distFromT, sr*5)
        # else:
        #     lowest_MVCost = 0
        #     mv_tree = None

        connection_cost = TotalTransformerCost + (sum(LVCostDict.values())) * LV
        # total_cost = TotalTransformerCost + (sum(LVCostDict.values())) * LV + lowest_MVCost
        print(i)
        print('Connection cost ', connection_cost)
        print('Revenue ', revenue)
        # print('Worse case total cost', total_cost)
        outFile.write("%i %f\n" % (i, sumLVCostAtEachStep[len(centers)]))  # need to clean up

        all_ready_checked.append([baseClusterID, mergingClusterID])
        # all_ready_checked.append([mergingClusterID,baseClusterID])
        if (connection_cost * fraction_recovered <= revenue):# & (lowest_MVCost < mv_wiggle_room):  # accept merge if cost is good
            minNodesByClusterID = copy.deepcopy(nodesByClusterID)
            # minCenters = copy.deepcopy(centers)
            # minLVCostDict = LVCostDict.copy()
            minClusterByNode = copy.deepcopy(clusterByNode)
            minTotalCost = connection_cost
            minLVCostSum = sumLVCostAtEachStep[len(centers)]
            # minTree = mv_tree
            minCenters_ST = copy.deepcopy(tx_locations)
            connectedNodesByClusterID = get_connected_nodes(tx_locations, nodesByClusterID)

        else:  # delete segment from evaluation
            centers = copy.deepcopy(tmp_centers)
            nodesByClusterID = copy.deepcopy(tmp_nodesByClusterID)
            clusterByNode = copy.deepcopy(tmp_clusterByNode)
            LVCostDict = copy.deepcopy(tmp_LVCostDict)
        segments = generateSegmentsDemand(centers, sr, all_ready_checked, roi_years, cost_per_kwh, TCost, LV,
                                          fraction_recovered)

        # generate segments for remaining nodes
        try:  # to check if there is a segment on the graph or there is only one cluster  # bir tane break eden varsa bile devamini check ediyor!!!!!
            # seems this looks for the shortest segment with the lv less that distFromT
            minSeg = min(segments, key=lambda obj: obj.getWeight())  # finds points within sr with max demand
            flag, tempCenterX, tempCenterY = maxTempInClusterDist(minSeg, clusterByNode, nodesByClusterID, distFromT)
            #            maxDist, maxDemand ,tempCenterX, tempCenterY = maxTempInClusterDist(minSeg, clusterByNode, nodesByClusterID,distFromT) # find the largest distance to centroid of minSeg
            if not flag:  # if largest distance > 500 meters

                segments.sort(key=lambda obj: obj.getWeight())  # sort by smallest to largest distances
                for seg in segments:
                    if (seg.getWeight() > distFromT * 2):  # distance greater than 1000 sk
                        break
                    else:  # if distance is okay check if there are 2 closer points
                        #                        maxDist, maxDemand, tempCenterX, tempCenterY = maxTempInClusterDist(seg, clusterByNode, nodesByClusterID,distFromT)
                        flag, tempCenterX, tempCenterY = maxTempInClusterDist(seg, clusterByNode, nodesByClusterID,
                                                                              distFromT)
                        if flag:
                            minSeg = seg  ## identifies a new minSeg to go to if there is still room to add to the LV
                            break  # finds the next minimum segment to go to
        except Exception as err:
            print(err)
            break

    # here redo prims with smallest network possible
    outFile.close()
    # import ipdb; ipdb.set_trace()
    optimized_mv_cost, minTree_ST, minCenters_ST_final, nodesByClusterID_ST_final = get_smallest_possible_MV_cost(minCenters_ST, MV, LV, minNodesByClusterID,
                                                                  minClusterByNode, distFromT, sr * 10000000000000)
    return minTotalCost + optimized_mv_cost, minTree_ST, minCenters_ST_final, connectedNodesByClusterID, minLVCostSum,nodesByClusterID_ST_final  # need a new dict, minNodesByClusterID_ST, minLVCostSum_ST


def addLVSeg(tree, centers, nodesByClusterID):  # single points line from the root
    SegID = 1000000

    for centerID in centers.keys():
        try:
            netID = tree.getNetID(centers[centerID])
        except:
            netID = 0
            tree._nodesByNetID[0] = []
            tree._network[netID] = []

        for node in nodesByClusterID[centerID]:
            length = ((node.getX() - centers[centerID].getX()) ** 2 +
                      (node.getY() - centers[centerID].getY()) ** 2) ** (.5)
            newSeg = network.Seg(SegID, node, centers[centerID], length)
            tree._netIDByNode[node] = netID
            tree._nodesByNetID[netID].append(node)
            tree._network[netID].append(newSeg)
    return tree


def writeLVDictToText(statsFile, Dict):
    'Writes LVCostDict to a text file for batchPrimsForTransformers.py.'
    outFile = open(statsFile, "w")
    for key in Dict.keys():
        LVCost = Dict[key] * 10
        outFile.write("%(key)i %(LVCost)f\n" % vars())
    outFile.close()
    return 0


def writeCenterSizeToText(statsFile, Dict):
    outFile = open(statsFile, "w")
    for key in Dict.keys():
        size = Dict[key].getWeight()
        outFile.write("%(size)i \n" % vars())
    outFile.close()
    return 0


def get_relevant_grids(txtpath, batch_num, total_num_batches):
    grid_files = []
    txt_files = os.listdir(txtpath)
    for txt in txt_files:
        with open(txtpath + txt) as src:
            content = src.readlines()
            for c in content:
                grid_files.append(c.strip())

    print("Number of grid files found", len(grid_files))
    step = int(np.ceil(len(grid_files) / total_num_batches))
    start = int(step * batch_num)
    stop = (start + step) - 1
    mygrids = grid_files[start:stop]
    return mygrids


def get_debug_grids(txtpath, batch_num, total_num_batches):
    grid_files = []
    txt_files = os.listdir(txtpath)
    for txt in txt_files:
        with open(txtpath + txt) as src:
            content = src.readlines()
            for c in content:
                grid_files.append(c.strip())
    print("Number of grid files found", len(grid_files))
    mygrids = grid_files[int(batch_num):int(batch_num + 1.0)]
    return mygrids


def get_debug_subgrids(txtpath, batch_num):
    # grid_files = []
    ward_files = os.listdir(txtpath)
    grid_files = sorted(ward_files)
    print("Number of grid files found", len(grid_files))
    mygrids = grid_files[int(batch_num):int(batch_num + 1.0)]
    mygrids = [os.path.join(txtpath, m) for m in mygrids]
    return mygrids


def get_relevant_wards(path, batch_num, total_num_batches):
    wards_file = os.listdir(path)
    wards = [os.path.join(path, c) for c in wards_file]
    print("Number of grid files found", len(wards))
    step = int(np.ceil(len(wards) / total_num_batches))
    start = int(step * batch_num)
    stop = (start + step) - 1
    mywards = wards[start:stop]
    return mywards


def check_start(start):
    val = start - np.floor(start)
    diffs = [val - 0, 1 - val]
    flag = np.argmin(diffs)
    if flag == 1:
        return int(np.ceil(start))
    else:
        return int(np.floor(start))


def get_scale_and_subgrids(grid, start, stop=None):
    # allowed_subgrid = np.arange(start,stop+1)
    subgrids = []
    for root, dirs, files in os.walk(grid):
        for name in files:
            if 'MV' not in name and 'FinalGrid' not in name:
                if 'scale' in root and '.shp' in name:
                    valid = int(name[:-4].split('_')[-1])
                    if valid == start:
                        subgrids.append(os.path.join(root, name))
                        return subgrids


def get_scale_and_subgrids_v2(grid, start, stop):
    # allowed_subgrid = np.arange(start,stop+1)
    subgrids = []
    for root, dirs, files in os.walk(grid):
        for name in files:
            if 'MV' not in name and 'FinalGrid' not in name:
                if 'scale' in root and '.shp' in name:
                    valid_stop = (name[:-4].split('/')[-1]).split('_')[-1]
                    valid_start = int((name[:-4].split('/')[-1]).split('_')[1])
                    if valid_start == start and valid_stop == stop:
                        subgrids.append(os.path.join(root, name))
                        return subgrids


def is_uncompleted(grid_file, subgrid_number):
    # files_in_subgrid = os.listdir(sub_grid)
    valid = True
    for root, dirs, files in os.walk(grid_file):
        for name in files:
            if 'scale' in root and 'modelOutput.txt' in name:
                cur_subgrid = int((name.split('_')[-1]).strip('modelOutput.txt'))
                if cur_subgrid == subgrid_number:
                    valid = False
    return valid


def is_notvisited(grid_file):
    valid = True
    for root, dirs, files in os.walk(grid_file):
        for name in files:
            if 'modelOutput.txt' in name:
                valid = False
    return valid


def get_search_radius(input_file, buffer_r=1000):
    bbox = fiona.open(input_file).bounds
    radius = int(np.ceil((((bbox[2] - bbox[0]) ** 2) + ((bbox[3] - bbox[1]) ** 2)) ** 0.5) + buffer_r)
    return radius


def convert_to_utm(latlon_shp, output_file):
    df = gpd.read_file(latlon_shp)
    df_utm = df.to_crs({'init': 'epsg:31028'})
    df_utm.to_file(output_file)


def main(cur_file):
    searchRadius = 1000  # meters
    # Cost parameters:
    MV = 25  # Cost of MV per meter
    LV = 10  # Cost of LV per meter
    TCost = 2000 / 8800.0  # Transformer Cost per watt
    distFromT = 500  # Dmax, direct distance from transformers
    maxLVLenghtInCluster = 500  # Lmax
    txCap = 8.8  # assuming 11kVA TX at a pf 0.8 to give 8.8 kW
    #    investment_per_hh = 300 # cost per HH
    investment = 600000
    roi_years = 5
    cost_per_kwh = 0.05  # although it is 0.20 a kwh recovery only about 5 cents
    fraction_covered_by_subsidy = 0
    fraction_to_be_recovered = 1 - fraction_covered_by_subsidy
    offgrid_cost = 500 ## assume offgrid cost
    # read shape file
    outputDir = cur_file[:-4]
    if os.path.isdir(outputDir):
        shutil.rmtree(outputDir)
    print(outputDir)
    logfilename = outputDir + 'modelStatus.txt'
    startTime = time.time()
    try:
        print("Generating Dictionaries")
        nodesByClusterID, clusterByNode, nodes, centers, LVCostDict, node_weights_output = generateDictsFromShp(
            cur_file,
            outputDir)

        node_weights_output.to_csv(os.path.join(outputDir, 'node_weights.csv'))
        print("Run function starts...")
        totalCost, tree, centers, connected_nodesByClusterID, LVCostSum, nodesByClusterID = run(centers, nodesByClusterID, clusterByNode,
                                                                    LVCostDict, searchRadius, MV, LV, TCost,
                                                                    distFromT, investment, cost_per_kwh, roi_years,
                                                                    fraction_to_be_recovered,
                                                                    maxLVLenghtInCluster, outputDir, logfilename)
        fileRW.genShapefile(tree, outputDir + ".prj", outputDir + os.sep + "MV.shp")

        # import ipdb;ipdb.set_trace()

        statsFile1 = outputDir + os.sep + "LVCostDict.txt"
        statsFile2 = outputDir + os.sep + "CenterSize.txt"
        writeLVDictToText(statsFile1, LVCostDict)
        writeCenterSizeToText(statsFile2, centers)
        MVLength = tree.getTotalEdgeWeight()
        MVCost = MVLength * MV
        numTransformer = len(centers)
        try:
            netID = tree.getNetID([*centers.values()][0])
        except:
            netID = 0
            tree._nodesByNetID[0] = []
            tree._network[netID] = []
        my_lv = 0
        connected_nodes = [len(nodesByClusterID[k]) for k in centers.keys()]
        connected_nodes = sum(connected_nodes)
        for ID in centers.keys():
            nodesByNodeID = {}
            segments, lvCost = CMST_dfs_OLD.CMST(nodesByClusterID[ID], maxLVLenghtInCluster, centers[ID])
            my_lv += lvCost
            for segment in segments.values():
                node1 = segment.getNode1()
                node2 = segment.getNode2()
                if node1.getID() not in nodesByNodeID.keys():
                    nodesByNodeID[node1.getID()] = node1
                if node2.getID() not in nodesByNodeID.keys():
                    nodesByNodeID[node2.getID()] = node2

            for node in nodesByNodeID.values():
                tree._netIDByNode[node] = netID
                tree._nodesByNetID[netID].append(node)

            for segment in segments.values():
                tree._network[netID].append(segment)

        fileRW.genShapefile(tree, outputDir + ".prj", outputDir + os.sep + "FinalGrid.shp")
        with open(outputDir + 'modelOutput.txt', 'w') as dst:
            dst.write("NumStructures:" + str(len(nodes)) + "\n")
            dst.write("NumConnectedStructures:" + str(connected_nodes) + "\n")
            dst.write("LVLength:" + str(my_lv) + "\n")
            dst.write("LVPerCustomer:" + str(float(my_lv) / connected_nodes) + "\n")
            dst.write("MVLength:" + str(MVLength) + "\n")
            dst.write("MVPerCustomer:" + str(MVLength / connected_nodes) + "\n")
            dst.write("Num Transformers:" + str(numTransformer) + "\n")
            dst.write("Customers Per Tx:" + str(connected_nodes / float(numTransformer)) + "\n")
            dst.write("Total LV Cost:" + str(my_lv * float(LV)) + "\n")
            dst.write("Total MV Cost:" + str(MVCost) + "\n")
            transformerCost = numTransformer * TCost
            dst.write("Transformer Cost:" + str(transformerCost) + "\n")
            total_cost = MVCost + my_lv * float(LV) + transformerCost
            dst.write("Total Grid Cost:" + str(total_cost) + "\n")
            dst.write("Offgrid Cost:" + str(offgrid_cost*(len(nodes)-connected_nodes)) + "\n")
            runningT = time.time() - startTime
            dst.write("Total Running Time:" + str(runningT) + "\n")
            runningT1 = time.time() - startTime
            dst.write("Final Running Time:" + str(runningT1))
    except Exception as error_with_grid:
        print("Error with file:", error_with_grid)


if __name__ == "__main__":
    #    my_file = '../../yuezi/senegal_gabar_site1_centroids_wArea_from_KMZ_JB_EPSG31028/senegal_gabar_site1_centroids_wArea_from_KMZ_JB_EPSG31028.shp'
    # my_file = '../../yuezi/discrete_SoyNorth_scale3_subgrid2_subgrid_2/subgrid_2.shp'
    my_file = '../../yuezi/mixed_demand_SoyNorth_scale3_subgrid2_subgrid_2/subgrid_2.shp'
    main(my_file)
