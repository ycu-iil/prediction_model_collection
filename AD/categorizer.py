import numpy as np
from rdkit import DataStructs
from sklearn.metrics.pairwise import cosine_similarity

def tanimoto_similarity_based(fp, ref_fps, conf):
    sim = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
    cutoff = conf["AD"]["cutoff"]
    Nmin = conf["AD"]["Nmin"]
    n = sum(x >= cutoff for x in sim)
    if n >= Nmin:
        return True
    else:
        return False

def euclidean_distance_based(data, refs, conf):
    sim = []
    for ref in refs:
        distance = np.linalg.norm(ref - data)
        sim.append(1 / (1 + distance))
    cutoff = conf["AD"]["cutoff"]
    Nmin = conf["AD"]["Nmin"]
    n = sum(x >= cutoff for x in sim)
    if n >= Nmin:
        return True
    else:
        return False

def cosine_similarity_based(data, refs, conf):
    data_2d = data.reshape(-1, len(data))   # Convert to a one-row 2-D array
    sim = cosine_similarity(data_2d, refs)
    cutoff = conf["AD"]["cutoff"]
    Nmin = conf["AD"]["Nmin"]
    n = sum(x >= cutoff for x in sim[0])   # Need [0] because `sim` is a one-row 2-D array
    if n >= Nmin:
        return True
    else:
        return False