from gadata import gadata
import numpy as np
import MDSplus as mds


def avg_ped(shot, tmin, tmax):
    conn = mds.Connection('localhost')
    def get_avg(tag):
        ped_obj = gadata(tag, shot, connection=conn)
        min_idx = np.where(ped_obj.xdata > tmin)[0].min()
        max_idx = np.where(ped_obj.xdata < tmax)[0].max()
        avg_ped = ped_obj.zdata[min_idx:max_idx].mean()
        return avg_ped

    avg_neped = get_avg('prmtan_neped')
    avg_teped = get_avg('prmtan_teped')
    avg_tribot = get_avg('tribot')
    avg_tritop = get_avg('tritop')

    print('ne ped: {:.3e}'.format(avg_neped))
    print('Te ped: {:.3f}'.format(avg_teped))
    print('Min. tri: {:.3f}'.format(min(avg_tribot, avg_tritop)))
