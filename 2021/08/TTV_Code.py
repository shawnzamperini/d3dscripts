import EFIT.load_gfile_d3d as loadg
import MDSplus as MDS
from scipy.io import readsav
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl


def main(shot, carbon_type='C2', time = 1900, server = 'localhost', tree='EFIT01', verbose=True):
    MDSplusCONN = MDS.Connection(server)
    parmDICT = loadg.read_g_file_mds(shot = shot, time=time, tree=tree, connection=MDSplusCONN, write2file=False)
    if carbon_type == C2:
        fromIDL_after_puff = readsav('emission_structure_out_cam0par_'+str(shot)+'_before_puff.sav')
        fromIDL_after_puff = readsav('emission_structure_out_cam0par_'+str(shot)+'_after_puff.sav')
        inverted_b4 = fromIDL_b4_puff('emission_structure')['inverted'][0]
        inverted_after = fromIDL_after_puff('emission_structure')['inverted'][0]
        if verbose:
            print(parmDICT)
            print(parmDICT.keys())
        inverted_b4_avg = np.mean()

def runthis():
    connection = MDS.Connection('localhost')
    tree = 'EFIT01'
    shot = 184272

    time_b4 = 2000
    time_after = 3900

    parmDICT_b4 = loadg.read_g_file_mds(shot = shot, time=time_b4, tree=tree, connection=connection, write2file=False)
    Rs_b4, Zs_b4 = np.meshgrid(parmDICT_b4['R'], parmDICT_b4['Z'])
    parmDICT_after = loadg.read_g_file_mds(shot = shot, time=time_after, tree=tree, connection=connection, write2file=False)
    Rs_after, Zs_after = np.meshgrid(parmDICT_after['R'], parmDICT_after['Z'])

    fromIDL_c2_b4 = readsav('emission_structure_out_cam0par_'+str(shot)+'_f120_f135.sav')
    fromIDL_c2_after = readsav('emission_structure_out_cam0par_'+str(shot)+'_f234_f246.sav')
    fromIDL_c3_b4 = readsav('emission_structure_out_cam0perp_'+str(shot)+'_f120_f135.sav')
    fromIDL_c3_after = readsav('emission_structure_out_cam0perp_'+str(shot)+'_f234_f246.sav')

    emission_c2_b4 = fromIDL_c2_b4('emission_structure')
    emission_c2_after = fromIDL_c2_after('emission_structure')
    emission_c3_b4 = fromIDL_c3_b4('emission_structure')
    emission_c3_after = fromIDL_c3_after('emission_structure')

    inverted_c2_b4 = emission_c2_b4['inverted'][0]
    inverted_c2_after = emission_c2_after['inverted'][0]
    inverted_c3_b4 = emission_c3_b4['inverted'][0]
    inverted_c3_after = emission_c3_after['inverted'][0]

    inverted_c2_b4_avg = (inverted_c2_b4[0]+inverted_c2_b4[1]+inverted_c2_b4[2]+inverted_c2_b4[3]+inverted_c2_b4[4]+
                            inverted_c2_b4[5]+inverted_c2_b4[6]+inverted_c2_b4[7]+inverted_c2_b4[8]+inverted_c2_b4[9]+
                            inverted_c2_b4[10]+inverted_c2_b4[11]+inverted_c2_b4[12]+inverted_c2_b4[13]+inverted_c2_b4[14])/len(inverted_c2_b4)
    inverted_c2_after_avg = (inverted_c2_after[0]+inverted_c2_after[1]+inverted_c2_after[2]+inverted_c2_after[3]+inverted_c2_after[4]+
                            inverted_c2_after[5]+inverted_c2_after[6]+inverted_c2_after[7]+inverted_c2_after[8]+inverted_c2_after[9]+
                            inverted_c2_after[10]+inverted_c2_after[11])/len(inverted_c2_after)
    inverted_c3_b4_avg = (inverted_c3_b4[0]+inverted_c3_b4[1]+inverted_c3_b4[2]+inverted_c3_b4[3]+inverted_c3_b4[4]+
                            inverted_c3_b4[5]+inverted_c3_b4[6]+inverted_c3_b4[7]+inverted_c3_b4[8]+inverted_c3_b4[9]+
                            inverted_c3_b4[10]+inverted_c3_b4[11]++inverted_c3_b4[12]+inverted_c3_b4[13]+inverted_c3_b4[14])/len(inverted_c3_b4)
    inverted_c3_after_avg = (inverted_c3_after[0]+inverted_c3_after[1]+inverted_c3_after[2]+inverted_c3_after[3]+inverted_c3_after[4]+
                            inverted_c3_after[5]+inverted_c3_after[6]+inverted_c3_after[7]+inverted_c3_after[8]+inverted_c3_after[9]+
                            inverted_c3_after[10]+inverted_c3_after[11])/len(inverted_c3_after)

    c2_diff = inverted_c2_after_avg - inverted_c2_b4_avg
    c3_diff = inverted_c3_after_avg - inverted_c3_after_avg

    im = []
    im.append(inverted_c2_b4_avg)
    im.append(inverted_c3_b4_avg)
    im.append(inverted_c2_after_avg)
    im.append(inverted_c3_after_avg)
    im.append(c2_diff)
    im.append(c3_diff)
    for i in range(6):
        print(i, '\n')
        print(im[i])

    title = []
    t1 = 'C2 before puff'
    t2 = 'C3 before puff'
    t3 = 'C2 after puff'
    t4 = 'C3 after puff'
    t5 = 'C2 diff'
    t6 = 'C3 diff'
    title.append(t1)
    title.append(t2)
    title.append(t3)
    title.append(t4)
    title.append(t5)
    title.append(t6)

    parm_wall1 = []
    parm_wall2 = []
    parm_lcfs1 = []
    parm_lcfs2 = []
    parm_psi = []
    R = []
    Z = []
    parm_wall1.append(parmDICT_b4['wall'][:, 0])
    parm_wall1.append(parmDICT_b4['wall'][:, 0])
    parm_wall1.append(parmDICT_after['wall'][:, 0])
    parm_wall1.append(parmDICT_after['wall'][:, 0])
    parm_wall1.append(parmDICT_after['wall'][:, 0]) #for_diff_plot (FDP)
    parm_wall1.append(parmDICT_after['wall'][:, 0]) #for_diff_plot

    parm_wall2.append(parmDICT_b4['wall'][:, 1])
    parm_wall2.append(parmDICT_b4['wall'][:, 1])
    parm_wall2.append(parmDICT_after['wall'][:, 1])
    parm_wall2.append(parmDICT_after['wall'][:, 1])
    parm_wall2.append(parmDICT_after['wall'][:, 1]) #FDP
    parm_wall2.append(parmDICT_after['wall'][:, 1]) #FDP

    R.append(Rs_b4)
    R.append(Rs_b4)
    R.append(Rs_after)
    R.append(Rs_after)
    R.append(Rs_after) #FDP
    R.append(Rs_after) #FDP

    Z.append(Zs_b4)
    Z.append(Zs_b4)
    Z.append(Zs_after)
    Z.append(Zs_after)
    Z.append(Zs_after) #FDP
    Z.append(Zs_after) #FDP

    parm_psi.append(parmDICT_b4['psiRZn'])
    parm_psi.append(parmDICT_b4['psiRZn'])
    parm_psi.append(parmDICT_after['psiRZn'])
    parm_psi.append(parmDICT_after['psiRZn'])
    parm_psi.append(parmDICT_after['psiRZn']) #FDP
    parm_psi.append(parmDICT_after['psiRZn']) #FDP

    parm_lcfs1.append(parmDICT_b4['lcfs'][:, 0])
    parm_lcfs1.append(parmDICT_b4['lcfs'][:, 0])
    parm_lcfs1.append(parmDICT_after['lcfs'][:, 0])
    parm_lcfs1.append(parmDICT_after['lcfs'][:, 0])
    parm_lcfs1.append(parmDICT_after['lcfs'][:, 0]) #FDP
    parm_lcfs1.append(parmDICT_after['lcfs'][:, 0]) #FDP

    parm_lcfs2.append(parmDICT_b4['lcfs'][:, 1])
    parm_lcfs2.append(parmDICT_b4['lcfs'][:, 1])
    parm_lcfs2.append(parmDICT_after['lcfs'][:, 1])
    parm_lcfs2.append(parmDICT_after['lcfs'][:, 1])
    parm_lcfs2.append(parmDICT_after['lcfs'][:, 1]) #FDP
    parm_lcfs2.append(parmDICT_after['lcfs'][:, 1]) #FDP

    radii_min = []
    radii_max = []
    elevation_min = []
    elevation_max = []

    radii_min.append(min(emission_c2_b4['radii'][0][0]))
    radii_max.append(max(emission_c2_b4['radii'][0][0]))
    radii_min.append(min(emission_c3_b4['radii'][0][0]))
    radii_max.append(max(emission_c3_b4['radii'][0][0]))
    radii_min.append(min(emission_c2_after['radii'][0][0]))
    radii_max.append(max(emission_c2_after['radii'][0][0]))
    radii_min.append(min(emission_c3_after['radii'][0][0]))
    radii_max.append(max(emission_c3_after['radii'][0][0]))
    radii_min.append(min(emission_c2_after['radii'][0][0])) #FDP
    radii_max.append(max(emission_c2_after['radii'][0][0])) #FDP
    radii_min.append(min(emission_c3_after['radii'][0][0])) #FDP
    radii_max.append(max(emission_c3_after['radii'][0][0])) #FDP


    elevation_max.append(max(emission_c2_b4['elevation'][0][0]))
    elevation_min.append(min(emission_c2_b4['elevation'][0][0]))
    elevation_max.append(max(emission_c3_b4['elevation'][0][0]))
    elevation_min.append(min(emission_c3_b4['elevation'][0][0]))
    elevation_max.append(max(emission_c2_after['elevation'][0][0]))
    elevation_min.append(min(emission_c2_after['elevation'][0][0]))
    elevation_max.append(max(emission_c3_after['elevation'][0][0]))
    elevation_min.append(min(emission_c3_after['elevation'][0][0]))
    elevation_max.append(max(emission_c2_after['elevation'][0][0])) #FDP
    elevation_min.append(min(emission_c2_after['elevation'][0][0])) #FDP
    elevation_max.append(max(emission_c3_after['elevation'][0][0])) #FDP
    elevation_min.append(min(emission_c3_after['elevation'][0][0])) #FDP

    fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True)
    i=0
    for ax in axes.flat:
        im1 = ax.imshow(im[i][:][:], cmap='magma', interpolation='quadric', extent=[radii_min[i], radii_max[i], elevation_max[i], elevation_min[i]])
        ax.set_title(title[i])
        ax.invert_yaxis()
        ax.plot(parm_wall1[i], parm_wall2[i], 'w--', linewidth=2)
        ax.axis([0.98, 1.6, 0.8, 1.4])
        ax.contour(R[i], Z[i], parm_psi[i], [0.9996], colors='white', linewidths=1)
        ax.plot(parm_lcfs1[i], parm_lcfs2[i], 'white', linewidth=1)
        ax.set_ylabel('Z[m]')
        ax.set_xlabel('R[m]')
        i+=1


    cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
    plt.colorbar(im1, cax=cax, **kw)
    #plt.set_ticks(np.arange(0, 1.1, 0.5))
    #plt.set_ticklabels(['low', 'medium', 'high'])
    plt.suptitle('Shot # '+str(shot))
    plt.show()

    #plt.figure()
    #plt.imshow(c2_diff, cmap='magma')
    #plt.title('C2_diff')

    #plt.figure()
    #plt.imshow(c3_diff, cmap='magma')
    #plt.title('C3_DIff')
    """
    #fig, axes = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True)
    #plotting C2_before
    plt.figure()
    ax1 = plt.subplot(231)
    im1 = ax1.imshow(inverted_c2_b4_avg, cmap='magma', interpolation='quadric', extent=[min(emission_c2_b4['radii'][0][0]), max(emission_c2_b4['radii'][0][0]), max(emission_c2_b4['elevation'][0][0]), min(emission_c2_b4['elevation'][0][0])])
    ax1.invert_yaxis()
    ax1.plot(parmDICT_b4['wall'][:, 0], parmDICT_b4['wall'][:, 1], 'w--', linewidth=2)
    ax1.axis([0.98, 1.6, 0.8, 1.4])
    ax1.contour(Rs_b4, Zs_b4, parmDICT_b4['psiRZn'], [0.9996], colors='white', linewidths=1)
    ax1.plot(parmDICT_b4['lcfs'][:, 0], parmDICT_b4['lcfs'][:, 1], 'white', linewidth=1)
    ax1.set_ylabel('Z[m]')
    ax1.text(1.05, 1.3, 'C2_before', c='white')
    plt.colorbar(im1)
    #plotting C2_after
    ax2 = plt.subplot(232)
    im2 = ax2.imshow(inverted_c2_after_avg, cmap='magma', interpolation='quadric', extent=[min(emission_c2_after['radii'][0][0]), max(emission_c2_after['radii'][0][0]), max(emission_c2_after['elevation'][0][0]), min(emission_c2_after['elevation'][0][0])])
    ax2.invert_yaxis()
    ax2.plot(parmDICT_after['wall'][:, 0], parmDICT_after['wall'][:, 1], 'w--', linewidth=2)
    ax2.axis([0.98, 1.6, 0.8, 1.4])
    ax2.contour(Rs_after, Zs_after, parmDICT_after['psiRZn'], [0.9996], colors='white', linewidths=1)
    ax2.plot(parmDICT_after['lcfs'][:, 0], parmDICT_after['lcfs'][:, 1], 'white', linewidth=1)
    ax2.set_ylabel('Z[m]')
    ax2.text(1.05, 1.3, 'C2_after', c='white')
    plt.colorbar(im2)
    #plotting C2_diff
    ax3 = plt.subplot(233)
    im3 = ax3.imshow(c2_diff, cmap='magma', interpolation='quadric', extent=[min(emission_c2_after['radii'][0][0]), max(emission_c2_after['radii'][0][0]), max(emission_c2_after['elevation'][0][0]), min(emission_c2_after['elevation'][0][0])])
    ax3.invert_yaxis()
    ax3.plot(parmDICT_after['wall'][:, 0], parmDICT_after['wall'][:, 1], 'w--', linewidth=2)
    ax3.axis([0.98, 1.6, 0.8, 1.4])
    ax3.contour(Rs_after, Zs_after, parmDICT_after['psiRZn'], [0.9996], colors='white', linewidths=1)
    ax3.plot(parmDICT_after['lcfs'][:, 0], parmDICT_after['lcfs'][:, 1], 'white', linewidth=1)
    ax3.set_ylabel('Z[m]')
    ax3.text(1.05, 1.3, 'C2_Diff', c='white')
    plt.colorbar(im3)

    #plotting C3_before
    ax4 = plt.subplot(234)
    im4 = ax4.imshow(inverted_c3_b4_avg, cmap='magma', interpolation='quadric', extent=[min(emission_c3_b4['radii'][0][0]), max(emission_c3_b4['radii'][0][0]), max(emission_c3_b4['elevation'][0][0]), min(emission_c3_b4['elevation'][0][0])])
    ax4.invert_yaxis()
    ax4.plot(parmDICT_b4['wall'][:, 0], parmDICT_b4['wall'][:, 1], 'w--', linewidth=2)
    ax4.axis([0.98, 1.6, 0.8, 1.4])
    ax4.contour(Rs_b4, Zs_b4, parmDICT_b4['psiRZn'], [0.9996], colors='white', linewidths=1)
    ax4.plot(parmDICT_b4['lcfs'][:, 0], parmDICT_b4['lcfs'][:, 1], 'white', linewidth=1)
    ax4.set_ylabel('Z[m]')
    ax4.set_xlabel('R[m]')
    ax4.text(1.05, 1.3, 'c3_before', c='white')
    plt.colorbar(im4)
    #plotting c3_after
    ax5 = plt.subplot(235)
    im5 = ax5.imshow(inverted_c3_after_avg, cmap='magma', interpolation='quadric', extent=[min(emission_c3_after['radii'][0][0]), max(emission_c3_after['radii'][0][0]), max(emission_c3_after['elevation'][0][0]), min(emission_c3_after['elevation'][0][0])])
    ax5.invert_yaxis()
    ax5.plot(parmDICT_after['wall'][:, 0], parmDICT_after['wall'][:, 1], 'w--', linewidth=2)
    ax5.axis([0.98, 1.6, 0.8, 1.4])
    ax5.contour(Rs_after, Zs_after, parmDICT_after['psiRZn'], [0.9996], colors='white', linewidths=1)
    ax5.plot(parmDICT_after['lcfs'][:, 0], parmDICT_after['lcfs'][:, 1], 'white', linewidth=1)
    ax5.set_ylabel('Z[m]')
    ax5.set_xlabel('R[m]')
    ax5.text(1.05, 1.3, 'c3_after', c='white')
    plt.colorbar(im5)
    #plotting c3_diff
    ax6 = plt.subplot(236)
    im6 = ax6.imshow(c3_diff, cmap='magma', interpolation='quadric', extent=[min(emission_c3_after['radii'][0][0]), max(emission_c3_after['radii'][0][0]), max(emission_c3_after['elevation'][0][0]), min(emission_c3_after['elevation'][0][0])])
    ax6.invert_yaxis()
    ax6.plot(parmDICT_after['wall'][:, 0], parmDICT_after['wall'][:, 1], 'w--', linewidth=2)
    ax6.axis([0.98, 1.6, 0.8, 1.4])
    ax6.contour(Rs_after, Zs_after, parmDICT_after['psiRZn'], [0.9996], colors='white', linewidths=1)
    ax6.plot(parmDICT_after['lcfs'][:, 0], parmDICT_after['lcfs'][:, 1], 'white', linewidth=1)
    ax6.set_ylabel('Z[m]')
    ax6.set_xlabel('R[m]')
    ax6.text(1.05, 1.3, 'c3_Diff', c='white')
    plt.colorbar(im6)

    plt.tight_layout()
    """
