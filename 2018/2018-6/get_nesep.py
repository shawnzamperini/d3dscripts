import get_ts as ts
import numpy  as np
import matplotlib.pyplot as plt



def get_nesep(shot, tmin=2500, tmax=4000, tstep=250):
    ts_dict = ts.run_script(shot, 'core', tmin=tmin, tmax=tmax, tstep=tstep)

    # Plot the data to make sure its not grabbing during an ELM.
    all_psins = ts_dict['psins']['all_psins']
    all_nes   = ts_dict['psins']['all_nes']
    all_tes   = ts_dict['psins']['all_Tes']
    plt.style.use('seaborn')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(all_psins, all_tes, '.')
    ax1.set_xlabel('Psin')
    ax1.set_ylabel('Te (eV)')
    ax1.set_title(str(shot))

    # Plot the avg values over it.
    ax1.plot(ts_dict['psins']['avg_psins'], ts_dict['psins']['avg_Tes'])
    plt.show()

    # Assuming it's a good fit, get the first point inside and outside sep.
    inside = np.min(np.where(ts_dict['psins']['avg_psins']<1.0))
    outside = inside - 1

    # Get avg ne between the two as ne_sep.
    ne_sep = np.mean([ts_dict['psins']['avg_nes'][inside], ts_dict['psins']['avg_nes'][outside]])
    Te_sep = np.mean([ts_dict['psins']['avg_Tes'][inside], ts_dict['psins']['avg_Tes'][outside]])

    print('ne sep = {:.4e}'.format(ne_sep))
    print('Te sep = {:.4e}'.format(Te_sep))

    return Te_sep
