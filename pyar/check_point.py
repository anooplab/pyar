import pickle


# import pandas as pd
# import numpy as np

def write_chcek_point(jobs_dict, location, logger):
    with open(location + '/jobs.pkl', 'wb') as f:
        pickle.dump(jobs_dict, f)
    return None


def read_check_point(location):
    flag = 0
    try:
        with open(location + '/jobs.pkl', 'rb') as f:
            x = pickle.load(f)
        for k in x.keys():
            if x[k] is not None:
                flag = 1
    except:
        flag = 0

    if flag == 1:
        return x
    if flag == 0:
        return None


def update_checkpoint(jobs_dict, to_remove, val, logger, location):
    if to_remove == 'gamma':
        jobs_dict.pop(val)
    if to_remove == 'ori':
        for i in jobs_dict.keys():
            for j in jobs_dict[i]:
                if j.name == val:
                    jobs_dict[i].remove(j)
    logger.debug("==========Updated Checkpoint: {}==========".format(val))
    write_chcek_point(jobs_dict, location, logger)
    return jobs_dict
