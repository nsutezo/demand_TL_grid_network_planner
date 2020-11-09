#!/usr/bin/env python
"""
A heuristic for Capacitated Minimum Spanning Tree (Esau-Williams Algorithm)
See page 8 at http://www.pitt.edu/~dtipper/2110/Slides6.pdf
Simone Fobi Akinola
"""

import numpy as np
from scipy.spatial import distance as dst

import network


def CMST_Caller(households, capacity, root):

    household_data = [(h.getID(), h.getX(), h.getY(), h.getWeight()) for h in households]
    household_data = np.array(household_data)
    root_data = [root.getID(), root.getX(), root.getY(), root.getWeight()]
    root_data = np.array(root_data)

    connections = CMST_only_rewrite(household_data, capacity, root_data)
    #print([(household_data[c[0],0],household_data[c[1],0]) for c in connections])

    top_parents = set(range(len(household_data))) - set([child for (child, parent) in connections])
    segments = [(int(household_data[child, 0]), int(household_data[parent, 0])) for (child, parent) in connections]
    segments += [(int(household_data[top_parent, 0]), root.getID()) for top_parent in top_parents]

    # total_LV_length = 0
    # for (child, parent) in connections:
    #     total_LV_length += np.linalg.norm(household_data[child, 1:3] - household_data[parent, 1:3])
    # for child in top_parents:
    #     total_LV_length += np.linalg.norm(household_data[child, 1:3] - root_data[1:3])
    # return segments, total_LV_length
    #import ipdb
    #ipdb.set_trace()
    connections += [(child, -root.getID() - 100) for child in top_parents]
    tree_segments = {}
    total_LV_length = 0
    #tmp_root =root
    tmp_root =network.Node((-root.getID() - 100), root.getX(), root.getY(), root.getWeight())

    for seg_id, (child, parent) in enumerate(connections):
        if parent < 0:
            length = np.linalg.norm(household_data[child, 1:3] - root_data[1:3])
            tree_segments[(int(household_data[child, 0]), parent)] = network.Seg(
                seg_id + int(1e7), households[child], tmp_root, length)
        else:
            length = np.linalg.norm(household_data[child, 1:3] - household_data[parent, 1:3])
            tree_segments[(int(household_data[child, 0]), int(household_data[parent, 0]))] = network.Seg(
                seg_id + int(1e7), households[child], households[parent], length)
        total_LV_length += length

    return tree_segments, total_LV_length

import copy
def CMST_only(household_data, capacity, root_data):
    pairwise_distance = dst.squareform(dst.pdist(household_data[:, 1:3]))
    np.fill_diagonal(pairwise_distance, np.infty)

    distance_to_root = np.linalg.norm(household_data[:, 1:3] - root_data[1:3], axis=1)


    group_assignment = {i: [i] for i in range(len(household_data))}
    distance_to_parents = np.zeros((len(household_data), 1))

    connections = []
    my_connections =[]
    done = False

    while not done:
        # compute tradeoff
        trade_off = np.min(pairwise_distance, axis=0) - distance_to_root
        # find the smallest valid trade_off, then merge and break
        for t_idx in np.argsort(trade_off):
            if trade_off[t_idx] < 0:
                closest_neighbor_idx = pairwise_distance[t_idx, :].argmin()  # parent node
                distance_to_furthest_point = distance_to_root[closest_neighbor_idx] + \
                                             pairwise_distance[closest_neighbor_idx, t_idx] + \
                                             np.max(distance_to_parents[group_assignment[t_idx]])


                if distance_to_furthest_point <= capacity:
                    # merge
                    my_connections.append((household_data[t_idx,0],household_data[closest_neighbor_idx,0]))
                    connections.append((t_idx, closest_neighbor_idx))
                    distance_to_parents[group_assignment[t_idx]] = distance_to_parents[group_assignment[t_idx]] + \
                                                                   pairwise_distance[closest_neighbor_idx, t_idx]
                    # distance_to_parents[t_idx] = pairwise_distance[closest_neighbor_idx, t_idx]

                    new_group = group_assignment[t_idx] + group_assignment[closest_neighbor_idx]

                    for node in new_group:
                        group_assignment[node] = list(set(new_group))
                    distance_to_root[group_assignment[t_idx]] = distance_to_root[closest_neighbor_idx]

                    #pairwise_distance[t_idx, :] = np.infty
                    #pairwise_distance[:, t_idx] = np.infty
                    pairwise_distance[t_idx, closest_neighbor_idx] = np.infty
                    pairwise_distance[closest_neighbor_idx, t_idx] = np.infty
                    break
                #else:
                #    pairwise_distance[t_idx, closest_neighbor_idx] = np.infty
                #    pairwise_distance[closest_neighbor_idx, t_idx] = np.infty
                #    break
            else:
                done = True
                break

    for c in connections:
        if c[0]==c[1]:
            connections.remove(c)
    return connections


def CMST_only_rewrite(household_data, capacity, root_data):
    pairwise_distance = dst.squareform(dst.pdist(household_data[:, 1:3]))
    org_pw = dst.squareform(dst.pdist(household_data[:, 1:3]))
    np.fill_diagonal(pairwise_distance, np.infty)

    distance_to_root = np.linalg.norm(household_data[:, 1:3] - root_data[1:3], axis=1)
    distance_to_parents = np.zeros((len(household_data), 1))

    children = {i: [i] for i in range(len(household_data))}

    connections = []
    my_connections =[]
    done = False
    #if household_data.shape[0]>7:
    #    import ipdb
    #    ipdb.set_trace()
    while not done:
        # compute tradeoff
        trade_off = np.min(pairwise_distance, axis=0) - distance_to_root
        # find the smallest valid trade_off, then merge and break
        #for t_idx in np.argsort(trade_off):
        t_idx = np.argmin(trade_off)
        if trade_off[t_idx] < 0:
            closest_neighbor_idx = pairwise_distance[t_idx, :].argmin()  # parent node

                # of the nodes connected to the parent, who is closest to current node _t_idx

                # get the children of closest_neighbor_idx
                # find the children of closest_neighbor_idx with shortest distance to cur_node
                # this could also return the current closest_neighbor_idx

            distance_with_children_of_closest_neighbor_idx = org_pw[t_idx, children[closest_neighbor_idx]]
            for k in np.argsort(distance_with_children_of_closest_neighbor_idx):
                best_child_closest_neighbor_idx = children[closest_neighbor_idx][k]
                # check capacity  if meets update else loop through the for loop
                # if best_child_closest_neighbor_idx == closest_neighbor_idx and capacity is violated, update  & break

                if best_child_closest_neighbor_idx != closest_neighbor_idx:
                    distance_to_furthest_point = distance_to_root[closest_neighbor_idx] + \
                                         np.max(distance_to_parents[best_child_closest_neighbor_idx]) +\
                                         org_pw[best_child_closest_neighbor_idx, t_idx] + \
                                         np.max(distance_to_parents[children[t_idx]])
                else:
                    distance_to_furthest_point = distance_to_root[closest_neighbor_idx] + \
                                            pairwise_distance[closest_neighbor_idx, t_idx] + \
                                            np.max(distance_to_parents[children[t_idx]])

                if distance_to_furthest_point <= capacity:
                    # merge
                    my_connections.append((household_data[t_idx, 0], household_data[best_child_closest_neighbor_idx, 0]))
                    #print(my_connections)
                    connections.append((t_idx, best_child_closest_neighbor_idx))
                    distance_to_parents[children[t_idx]] = distance_to_parents[children[t_idx]] + \
                                                           org_pw[best_child_closest_neighbor_idx, t_idx] + \
                                                           np.max(distance_to_parents[best_child_closest_neighbor_idx])

                    new_group = children[t_idx] + children[closest_neighbor_idx]

                    if closest_neighbor_idx != best_child_closest_neighbor_idx:
                        new_sub_group = children[t_idx] + children[best_child_closest_neighbor_idx]
                        children[best_child_closest_neighbor_idx] = list(set(new_sub_group))

                    children[closest_neighbor_idx] = list(set(new_group))

                    #for node in new_group:
                    #    children[node] = list(set(new_group))
                    distance_to_root[children[t_idx]] = distance_to_root[closest_neighbor_idx]
                    pairwise_distance[t_idx,:] = np.infty
                    pairwise_distance[:, t_idx] = np.infty
                    #pairwise_distance[t_idx, best_child_closest_neighbor_idx] = np.infty
                    #pairwise_distance[best_child_closest_neighbor_idx, t_idx] = np.infty

                    break
        else:
            done = True
    return connections



    # while min_trade_off < 0:
    #     # compute tradeoff
    #     t_idx = np.argmin(trade_off)
    #     closest_neighbor_idx = pairwise_distance[t_idx, :].argmin()  # parent node
    #     distance_to_furthest_point = distance_to_root[closest_neighbor_idx] + \
    #                                 pairwise_distance[closest_neighbor_idx, t_idx] + \
    #                                 np.max(distance_to_parents[children[t_idx]])
    #
    #     t_idx_parents = parents[t_idx]
    #     # suggested parent - current parent of t_idx
    #     for t_idx_parent in np.argsort(t_idx_parents):
    #         compare_parents_dist = org_pw[t_idx,closest_neighbor_idx]-org_pw[t_idx,t_idx_parent]
    #
    #     #t_idx_parent = t_idx_parents[np.argmin(org_pw[t_idx,closest_neighbor_idx]-org_pw[t_idx,t_idx_parents])]
    #     # if closer to new parent compared to current parent
    #     #if compare_parents <0 :
    #         #update parent child association and delete from connections
    #         if compare_parents_dist < 0 and distance_to_furthest_point <= capacity:
    #         #if closest_neighbor_idx in connected_j and t_idx in connected_i:
    #             connections.append((t_idx_parent, closest_neighbor_idx))
    #             print(connections)
    #
    #             new_group = list(set(children[closest_neighbor_idx] + children[t_idx]))
    #             children[closest_neighbor_idx] = new_group
    #             parents[t_idx] = list(set(closest_neighbor_idx + parents[t_idx]))
    #
    #             for k in children[t_idx]:
    #                 parents[k] = closest_neighbor_idx
    #
    #             distance_to_root[new_group] = distance_to_root[closest_neighbor_idx]
    #
    #             # distance_to_parents[group_assignment[t_idx]] = distance_to_parents[group_assignment[t_idx]] + \
    #             #                                                        pairwise_distance[closest_neighbor_idx, t_idx]
    #
    #             distance_to_parents[children[t_idx]] = distance_to_parents[children[t_idx]] - new_pw[t_idx,t_idx_parent]+ \
    #                                                        new_pw[closest_neighbor_idx, t_idx]
    #             break
    #
    #         #pairwise_distance[t_idx, :] = np.infty
    #         #pairwise_distance[:, t_idx] = np.infty
    #
    #
    #     pairwise_distance[t_idx, closest_neighbor_idx] = np.infty
    #     pairwise_distance[closest_neighbor_idx, t_idx] = np.infty
    #
    #     trade_off = np.min(pairwise_distance, axis=0) - distance_to_root
    #     min_trade_off = trade_off[np.argmin(trade_off)]
    #
    # for c in connections:
    #     if c[0]==c[1]:
    #         connections.remove(c)


