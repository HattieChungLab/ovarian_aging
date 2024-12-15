# All the methods needed to calculate the scaled distance within follicle and corpus luteum objects
# by Sophia Szady
def distance(p1, p2):
    """ Calculates the distance between 2 points in 2d space
    Notes: This method assumes p1 and p2 are in the format [x,y]
    Args: 
    p1: point 1
    p2: point 2 
    Returns:
    dist: the linear distance between p1 and p2
    """
    dist = math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)
    return dist

def centroid(arr):
    """ Calculates the centroid of a given object
    Notes: This method assumes that arr[:, 0] contains the x coordinates
    and arr[:, 1] contains the y coordinates
    Args:
    arr: an array of all the x and y coordinates of a given object
    Returns:
    x_coord: the x coordinate of the centroid
    y_coord: the y coordinate of the centroid
    """
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    x_coord = sum_x/length
    y_coord = sum_y/length
    return x_coord, y_coord

def farthest_point(points):
    """ Calculates the farthest point from the centroid
    Notes: see centroid and distance notes for data assumptions
    Args:
    points: an array of the x and y coordinates of a given object
    Return:
    cent: the centroid of the object
    farthest: the farthest point from the centroid in the object
    """
    cent = centroid(points)
    farthest = None
    max_distance = -1
    for point in points:
        dist = distance(cent, point)
        if dist > max_distance:
            max_distance = dist
            farthest = point
    return cent, farthest

def get_object_centroid_edge(slidex, object_id):
    """ Calculates the centroid and radius of an object 
    Notes: see centroid and distance notes for data assumptions,
    the method assumes the object is round/circular
    Args:
    slidex: an anndata object containing spatial data
    object_id: the id of the object of interest
    Return:
    obj_centroid: the centroid of the object
    obj_radius: the radius of the object
    """
    obj_coords = slidex[slidex.obs['segment_id']==object_id].obsm['spatial']
    obj_centroid, obj_edge = farthest_point(obj_coords)
    obj_radius = distance(obj_centroid, obj_edge)
    return obj_centroid, obj_radius

def calc_immune_programs_radial_distance(ad_x, ids):
    """ Calculates the scaled distance within each object
    Notes: see centroid and distance notes for data assumptions,
    the method assumes each object is round/circular
    Args:
    ad_x: an anndata object containing spatial data
    ids: a list of the ids of the objects of interest
    Return:
    ad_x: a modified version of the original anndata with the new obs 'scaled_dist'
    """
    for foll in ids: 
        idx_segment = ad_x.obs['segment_id']==foll
        obj_cent, obj_edge = get_object_centroid_edge(ad_x, foll)
        obj_coords = ad_x[idx_segment].obsm['spatial']
        # if the object is more than 1 point
        if len(obj_coords) > 1:
            # scaling the distance by the size of the object
            scaled_dist = [distance(obj_coords[i], obj_cent)/obj_edge for i in range(len(obj_coords))]
            ad_x.obs.loc[idx_segment, 'scaled_dist'] = scaled_dist
    return ad_x
